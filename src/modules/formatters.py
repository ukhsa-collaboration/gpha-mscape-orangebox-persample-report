import pandas as pd
import numpy as np
import logging

logger = logging.getLogger(__name__)

# add in the metadata fields
METADATA_FIELDS = {
    "climb_id": "CLIMB ID",
    "site": "Site", 
    "published_date": "Published date", 
    "collection_date": "Collection date", 
    "sample_source": "Sample source", 
    "sample_type": "Sample type"
}

HCID_COUNTS_COLUMNS = {"name", "taxon_id"}

SYLPH_COLUMNS = {
    "taxon_id" : "Taxon ID",
    "human_readable": "Species", 
}

SYLPH_TAXONOMY_COLUMNS = {
    "human_readable_species" : "Species",
    "gtdb_taxon_string" : "Taxonomy", 
    "gtdb_assembly_id" : "Assembly ID",
    "taxonomic_abundance" : "Abundance"
}

KRAKEN_COLUMNS = {
    "taxon_id": "Taxon ID", 
    "human_readable": "Species"
}

class FormattedSection:
    def __init__(self, status: str, table: pd.DataFrame | None = None, message: str = ""):
        self.status = status
        if table is None:
            table = pd.DataFrame()  # a new empty table for each section
        self.table = table
        self.message = message

def get_single_row(data: pd.DataFrame | None, section: str) -> pd.Series | None:
    """
    Return the only row of a per-sample DataFrame, None if there is no data, or raise if there are several lines.
    Enforces the assumption that each section of the report is per-sample, so there should be only one row of data.
    """
    if data is None or data.empty:
        return None
    if len(data) != 1:
        raise ValueError(f"{section}: expected data for one sample, got {len(data)} rows")
    return data.iloc[0]

def is_blank(value) -> bool:
    if value is None:
        return True
    if isinstance(value, str):
        return value.strip() == ""
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):  # lists and other containers
        return False

def format_date(value) -> str:
    """Format a date as DD/MM/YYYY. Raises ValueError if it cannot be parsed."""
    return pd.to_datetime(value, errors="raise").strftime("%d/%m/%Y")

def format_hcid(counts: pd.DataFrame | None, warnings:list[dict]) -> FormattedSection:
    """
    Input func deets here pls

    1. If any warning.json present -> flagged, list taxa and report. 
    2. Else hcid.counts.csv found -> screened and not flagged
    3. Else (nothing found) - Error or HCID check not run, reusult not available. 

    Args: 
        counts: the hcid.counts.csv table, or none if not found 
        warnings: the contents of each warning json
    """
    if counts is not None and not counts.empty:
        missing_columns = HCID_COUNTS_COLUMNS - set(counts.columns)
        if missing_columns:
            raise ValueError(f"Missing expected columns in hcid.counts.csv: {missing_columns}")

    # 1. Flaggeed
    if warnings:
        return format_hcid_flagged(counts, warnings)

    # 2. Screened and not flagged - phew
    if counts is not None and not counts.empty:
        logger.info("HCID check completed: sample screened and not flagged.")
        return format_hcid_screened(counts)

    # 3. Not run or an error raised
    logger.warning("HCID check not run or no results available.")
    return FormattedSection(status=STATUS_MISSING, message="HCID screening result not available for this sample.")

def format_hcid_flagged(counts: pd.DataFrame | None, warnings: list[dict]) -> FormattedSection:
    """
    Build the flagged message and table for the HCID section. 
    Look up taxon counts in the counts table if available, otherwise just report the taxon_id.
    """
    names = {} # taxon_id -> HCID name; such as 33278 = AI-A infection H7N9. 

    # find taxon and name for each warning and add to the table. 
    if counts is not None and not counts.empty: 
        names = dict(zip(counts["taxon_id"], counts["name"]))
    else:
        logger.error("HCID counts table not available, will report taxon_id only for flagged taxa.")

    rows = []
    for warning in warnings:
        if "taxon_id" not in warning:
            raise ValueError(f"Warning JSON missing 'taxon_id': {warning}")
        taxon_id = warning["taxon_id"]
        name = names.get(taxon_id)

        if name is None:
            # register a positive result, we can't find the parsewd info
            logger.error(f"HCID counts table does not contain taxon_id {taxon_id} from warning JSON.")
            name = f"Unknown taxon_id {taxon_id}"
        rows.append({"HCID taxon": taxon_id, "name": name})

    table = pd.DataFrame(rows)
    message = f"HCID check completed: sample flagged for {len(warnings)} taxon(s)."
    return FormattedSection(status=STATUS_OK, table=table, message=message)

def format_metadata(data: pd.DataFrame) -> FormattedSection:
        """
        Function to format the metadata of the run. Included are the following fields:
        •	Climb-id
        •	Site
        •	Published date
        •	Collection date
        •	Sample source
        •	Sample type


        Args: 
            data (pd.DataFrame): dataframe containing the flattened json metadata.

        Returns: 
            formatted section containing all of the sample metadata. 
        """
    if data is not None and not data.empty:
        row = get_single_row(data, "metadata")
    else:
        logger.warning("metadata is not present in the sample output")
        return FormattedSection(status = STATUS_MISSING, message="Sample metadata not found for this sample")

    rows = []
    # retrieve the metadata fields
    for field_name, label in METADATA_FIELDS.items():
        value = row.get(field_name)
        if is_blank(value):
            # record blank in logger
            logger.info(f"Metadata field {field_name} is missing or blank.")
            display = NOT_RECORDED
        elif field_name in METADATA_DATE_FIELDS:
            display = format_date(value)
        else:
            display = str(value)

    return FormattedSection(status=STATUS_OK, table=pd.DataFrame(rows, columns=["Field", "Value"]))
    
def format_kraken_genera(data: pd.DataFrame) -> FormattedSection:
    """
    Format the Kraken genera output for HTML and plotly output.

    Args:
        data (pd.DataFrame): pandas dataframe of kraken outputs
    """
    if data is None or data.empty:
        logger.warning(f"Kraken genus results not available")
        return FormattedSection(
            status = STATUS_MISSING, 
            message = "Kraken genus results not available for this sample."
        )
    
    missing_columns = KRAKEN_GENRA_COLUMNS - set(data.columns)
    if missing_columns:
        raise ValueError(f"Expected Kraken columns are missing")

    genera = data.loc[data["rank"] == "G", ["name", "taxon_id", "count"]].copy()
    genera["percentage"] = pd.to_numeric(genera["percentage"],errors="raise")
    top_genera = genera.sort_values(by="percentage", ascending=False).head(10)

    table = top_genera[["taxon_id", "human_readable", "percentage", "count_descendants", "profile"]]
    table.rename(columns={"taxon_id": "Taxon ID", "human_readable": "Name", "percentage": "Percentage", "count_descendants": "Descendant count", "profile": "Profile"}, inplace=True)

    message = f"Kraken genus results available for {len(genera)} genera, top 10 displayed."
    return FormattedSection(status=STATUS_OK, table=table, message=message)

def format_kraken_high_confidence_species(data: pd.DataFrame) -> FormattedSection:
    """
    Format the Kraken high confidence species output for HTML and plotly output.

    Args:
        data (pd.DataFrame): pandas dataframe of kraken outputs
    """
    if data is None or data.empty:
        logger.warning(f"Kraken high confidence species results not available")
        return FormattedSection(
            status = STATUS_MISSING, 
            message = "Kraken high confidence species results not available for this sample."
        )
    
    missing_columns = KRAKEN_COLUMNS - set(data.columns)
    if missing_columns:
        raise ValueError(f"Expected Kraken columns are missing")

    species = data.loc[data["rank"] == "S", ["name", "taxon_id", "count"]].copy()
    species["percentage"] = pd.to_numeric(species["percentage"],errors="raise")
    top_species = species.sort_values(by="percentage", ascending=False).head(10)

    table = top_species[["taxon_id", "human_readable", "percentage", "count_descendants", "profile"]]
    table.rename(columns={"taxon_id": "Taxon ID", "human_readable": "Name", "percentage": "Percentage", "count_descendants": "Descendant count", "profile": "Profile"}, inplace=True)

    message = f"Kraken high confidence species results available for {len(species)} species, top 10 displayed."
    return FormattedSection(status=STATUS_OK, table=table, message=message)

def format_sylph_results(data: pd.DataFrame) -> FormattedSelection:
    """
    Format the sylph output and get it ready for HTML and plotly display. 

    Args:
        Data (pd.Dataframe): pandas dataframe containing the sylph output.
    """
    if data is None or data.empty:
        logger.warning(f"High confidence Sylph results not available.")
        return FormattedSection(
            status = STATUS_MISSING,
            message = "Sylph high confidence results not available for this sample."
        )
    missing_columns = SYLPH_COLUMNS - set(data.columns)
    if missing_columns:
        raise ValueError(f"Columns missing")


    sylph_table = data.copy()
    # Fill in the blank human readable species columns with the human readable data where the human readable species value is NaN 
    sylph_table["species"] = sylph_table["human_readable"].fillna(sylph_table["human_readable"])

    sylph_table = sylph_table.sort_values(["effective_coverage", "adjusted_ani"], ascending = False)

    # isolate the high confidence species from sylph
    #TODO: check with Rachael is this is an ok way to do this/if anything should be returned for low confidence. 
    confidence = sylph_table["sylph_confidence"].astype(str).str.strip().str.lower()
    is_high = confidence == "high"
    n_high = int(sum(is_high))
    logger.info(f"Sylph results: {n_high} high confidence species found.")

    sylph_table = (sylph_table[list(SLYLPH_DISPLAY_COLUMNS)].rename(columns=SLYLPH_DISPLAY_COLUMNS).round(2).reset_index(drop=True))

    return FormattedSection(status=STATUS_OK, table=sylph_table, message=f"Sylph results available for {len(sylph_table)} species, {n_high} high confidence species found.")

SYLPH_TAXON = {"d" : "Domain", "p" : "Phylum", "c" : "Class", "o" : "Order", "f" : "Family", "g" : "Genus", "s" : "Species"}
RANK_COLUMNS = list(SYLPH_TAXON.values())

def split_sylph_lineage(lineages: pd.Series) -> pd.DataFrame:
    """
    Splits the gtdb_taxon_string column from sylph outputs to give taxonomic string. 

    """
    # check for nan and balnk cells
    blank = ""
    if lineages.isna().any() or (lineages.str.strip() == blank).any():
        logger.warning(f"na values in slyph taxonomy strings.")
        raise ValueError("Some sylph rows lack a taxonomy string")
    
    for lineage in lineages:
        row = {rank: f"Unclassified {rank.lower()}" for rank in RANK_COLUMNS}
        for part in lineage.split(";"):
            prefix, _, name = part.partition("__")
            if prefix in SYLPH_TAXON:
                row[SYLPH_TAXON[prefix]] = name.strip()
        rows.append(row)
    return pd.DataFrame(rows, columns=RANK_COLUMNS)

def format_taxonomy_sylph(data: pd.DataFrame) -> FormattedSelection:
    """
    This function will take the taxonomy string from the output sylph csv and return counts of each level.
    eg. d__Bacteria;p__Bacillota;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s__Faecalibacterium taiwanense
    - Domain = Bacteria
    - Phylum = Bacillota
    - Class = Clostridia
    - Order = Oscillospirales
    - Family = Ruminococcaceae
    - Genus = Faecalibacterium
    - Species = Faecalibacterium taiwanense
    # kingdom seems to be ignored

    Args: 
        data (pd.DataFrame): the sylph output csv
    
    Returns:
        counts (FormattedSelection): the output counts for each domain.
    """

    if data is None or data.empty:
        logger.warning(f"Sylph taxonomies unavailable from {data}.")
        return FormattedSection(status = STATUS_MISSING, message = "Slyph results unavailable for the selection.")

    missing_columns = SYLPH_TAXONOMY_COLUMNS - set(data.columns)
    if missing_columns:
        raise ValueError(f"key information not present")

    sylph_taxonomy = data.copy()
    sylph_taxonomy = split_sylph_lineage(sylph_taxonomy["gtdb_taxon_string"])
    # raise over coerce as we want nice and correct plotly sunbursts. 
    sylph_taxonomy["Taxonomic abundance (%)"] = pd.to_numeric(data["taxonomic_abundance"], errors="raise")

    logger.info("sylph taxonomy: %d genomes in %d phyla", len(sylph_taxonomy), sylph_taxonomy["Phylum"].nunique())
    return FormattedSection(status=STATUS_OK, table=sylph_taxonomy)



def summarize_sylph_taxonomy(data: pd.DataFrame | None) -> FormattedSection:
    """
    Split the Sylph lineages into one column per rank, ready for plotly sunbursts.

    Returns one row per genome (not counts): Domain ... Species and taxonomic abundance.
    Plotly adds the rows up, by abundance (values="Taxonomic abundance (%)") or by count (no values).

    Args:
        data (pd.DataFrame | None): the ClasPar Sylph table, or None if the file was not found.

    Returns:
        FormattedSection: the per-genome table, or a missing / no results message.
    """
    if data is None:
        logger.warning("sylph taxonomy: no ClasPar Sylph file, result unknown")
        return FormattedSection(status=STATUS_MISSING, message="Sylph results not available for this sample.")

    if data.empty:
        logger.info("sylph taxonomy: ClasPar Sylph file is empty")
        return FormattedSection(status=STATUS_NO_RESULTS, message="Sylph did not detect any bacteria or archaea.")

    missing_columns = SYLPH_TAXONOMY_COLUMNS - set(data.columns)
    if missing_columns:
        raise ValueError(f"Sylph data is missing columns: {sorted(missing_columns)}")

    # Number the rows 0, 1, 2 ... so they match the rows split_sylph_lineage returns.
    data = data.reset_index(drop=True)

    sylph_taxonomy = split_sylph_lineage(data["gtdb_taxon_string"])
    # raise over coerce: a bad abundance should stop the run, not quietly drop out of the sunburst.
    sylph_taxonomy["Taxonomic abundance (%)"] = pd.to_numeric(data["taxonomic_abundance"], errors="raise")

    logger.info("sylph taxonomy: %d genomes in %d phyla", len(sylph_taxonomy), sylph_taxonomy["Phylum"].nunique())
    return FormattedSection(status=STATUS_OK, table=sylph_taxonomy)


    

#TODO: add parser for the kraken sunburst
#TODO: add viral aligner formatter
















