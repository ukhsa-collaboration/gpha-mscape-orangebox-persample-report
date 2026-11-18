#!/usr/bin/env python3

import json
import pandas as pd
from pathlib import Path
import logging


def load_json(path: Path) -> dict:
    """
    Load a json file and return the contents as a dictionary. 
    
    Args: 
        path (Path): Path to the JSON file.

    Returns:
        dict: The contents of the JSON file as a dictionary.
    """
    try:
        with open(path, 'r') as f:
            data = json.load(f)
    except Exception as e:
        logging.error(f"Error loading JSON file from {path}: {e}")
        raise
    return data

def normalize_json(data: dict) -> pd.DataFrame:
    """
    Function takes dict data and normalizes it into a pandas DataFrame. 
    It also normalizes any nested fields specified in the methods list.

    Args:
        data (dict): The input data to be normalized.
    
    Returns:
        pd.DataFrame: A DataFrame containing the normalized data.
    """
    data = data.copy()

    for field in ("methods", "result_metrics"):
        if isinstance(data.get(field), str):
            data[field] = json.loads(data[field])

    return pd.json_normalize(data, sep="_")


def load_hcid_outputs(folder: Path) -> tuple[pd.DataFrame | None, list[dict]]:
    """
    Load the outputs of Scylla's check_hcid_status module for one sample.

    The folder contains:
        hcid.counts.csv (older samples: {CLIMB-ID}.hcid.counts.csv) - every screened sample has one
        {CLIMB-ID}.{taxid}.warning.json - one per flagged HCID taxon, flagged samples only

    Args:
        folder (Path): The sample's HCID output folder.

    Returns:
        The counts table (None if no counts file was found) and a list of warning dictionaries
        (empty if the sample did not flag).
    """
    folder = Path(folder)
    if not folder.is_dir():
        logging.warning(f"HCID folder not found: {folder}")
        return None, []

    counts_files = sorted(folder.glob("*hcid.counts.csv"))
    if not counts_files:
        logging.warning(f"No hcid.counts.csv found in {folder}")
        counts = None
    else:
        if len(counts_files) > 1:
            raise ValueError(f"More than one hcid.counts.csv found in {folder}: {counts_files}")
        counts = pd.read_csv(counts_files[0])
        logging.info(f"Loaded HCID counts from {counts_files[0]}")

    warning_files = sorted(folder.glob("*.warning.json"))
    warnings = [load_json(path) for path in warning_files]
    logging.info(f"Found {len(warnings)} HCID warning file(s) in {folder}")

    return counts, warnings
 
    

