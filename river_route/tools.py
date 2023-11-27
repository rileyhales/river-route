import pandas as pd


def configs_from_rapid(riv_bas_id: str,
                       k: str,
                       x: str,
                       rapid_connect: str,
                       out_params: str,
                       out_connectivity: str, ) -> None:
    """
    Generate river-route configuration files from input files for RAPID

    Args:
        riv_bas_id: Path to riv_bas_id CSV file
        k: Path to k CSV file
        x: Path to x CSV file
        rapid_connect: Path to rapid_connect CSV file
        out_params: Path to output routing parameters parquet file
        out_connectivity: Path to output connectivity parquet file

    Returns:
        None
    """
    pd.concat([
        pd.read_csv(riv_bas_id, header=None, names=['rivid']),
        pd.read_csv(x, header=None, names=['x']),
        pd.read_csv(k, header=None, names=['k']),
    ], axis=1).to_parquet(out_params)
    (
        pd
        .read_csv(rapid_connect, header=None)
        .iloc[:, :2]
        .rename(columns={0: 'rivid', 1: 'downstream_rivid'})
        .to_parquet(out_connectivity)
    )
    return
