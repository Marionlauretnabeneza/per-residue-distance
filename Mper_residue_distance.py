# Plotting and Saving to Excel
    
import json
from pathlib import PosixPath
import gemmi
from gemmi import cif, make_structure_from_block
import numpy as np
import pandas as pd
import logging
import argparse
from matplotlib import pyplot as plt


def calculate_static_rmsd(
    polymer1: gemmi.PolymerType.PeptideL,
    polymer2: gemmi.PolymerType.PeptideL,
    default_ptype: gemmi.PolymerType = gemmi.PolymerType.PeptideL,
) -> float:
    """
    Calculate the RMSD between two structures, without applying any transformation or
    superposition. This should have either been done already or intentionally skipped.
    The two structures must be of the same polymer type (left-handed protein by
    default).

    :param polymer1: First structure as Gemmi Polymer object
    :type polymer1: gemmi.PolymerType.PeptideL
    :param polymer2: Second structure as Gemmi Polymer object
    :type polymer2: gemmi.PolymerType.PeptideL
    :raises ValueError: Not protein structures
    :return: RMSD between the two structures
    :rtype: float
    """

    if polymer1.check_polymer_type() == polymer2.check_polymer_type() == default_ptype:
        rmsds = gemmi.calculate_current_rmsd(
            polymer1, polymer2, default_ptype, gemmi.SupSelect.CaP  # TODO Needed?
        )

        return rmsds.rmsd

    else:
        raise ValueError("Polymer types not the same or supported")




def load_mmcif(
    path: PosixPath | str, return_model: bool = False
) -> gemmi.cif.Block | gemmi.Structure:
    """
    Returns mmCIF Gemmi object given path to file.
    Modified from https://github.com/PDBeurope/protein-cluster-conformers

    :param path: Load (updated) mmCIF file into memory.
    :type path: pathlib.Path | str
    :param return_model: Return Gemmi Structure object instead of Block object
    :type return_model: bool
    :raises ValueError: File parsed was not mmCIF or gzip mmCIF
    :return: Contents of mmCIF as Gemmi Block object
    :rtype: gemmi.cif.Block | gemmi.Structure
    """
    path = str(path)

    if path[-6:] == "cif.gz":
        cif_block = cif.read(path).sole_block()

    elif path[-3:] == "cif":
        cif_block = cif.read_file(path).sole_block()

    else:
        raise ValueError("File parsed was not mmCIF or gzip mmCIF")

    # Convert to Gemmi Structure object if required
    return make_structure_from_block(cif_block)[0] if return_model else cif_block




def load_RT_matrices(path: PosixPath | str) -> dict:
    """
    Load a JSON file containing a set of RT matrices. The matrices included in the file
    should be 4x4 affine matrices, with the 3x3 rotation matrix in the upper left corner
    and the 3x1 translation vector in the right column.

    :param path: Path to JSON file
    :type path: PosixPath | str
    :return: Dictionary of RT matrices. Keyed by PDB and chain (auth_asym_id) ID.
    :rtype: dict
    """

    with open(path, "r") as f:
        rt_matrices = json.load(f)

    return rt_matrices





def affine_to_transformation(affine_matx: "list[float]") -> gemmi.Transform:
    """
    Convert an affine matrix to a Gemmi transformation object.

    :param affine_matx: 4x4 affine matrix. Includes the 3x3 rotation matrix in the upper
        left corner and the 3x1 translation vector in the right column.
    :type affine_matx: list[float]
    :return: Gemmi transformation object
    :rtype: gemmi.Transform
    """

    # Parse affine matrix into 3x3 rotation matrix
    m = [
        [affine_matx[0][0], affine_matx[0][1], affine_matx[0][2]],
        [affine_matx[1][0], affine_matx[1][1], affine_matx[1][2]],
        [affine_matx[2][0], affine_matx[2][1], affine_matx[2][2]],
    ]
    rotation_matrix = gemmi.Mat33(m)

    # Parse affine matrix into 3x1 translation vector
    translation_vector = gemmi.Vec3(
        affine_matx[0][3], affine_matx[1][3], affine_matx[2][3]
    )

    return gemmi.Transform(rotation_matrix, translation_vector)



def block_to_df(block: gemmi.cif.Block) -> pd.DataFrame:
    loop_names = [
        "group_PDB",
        "label_asym_id",
        "auth_asym_id",
        "pdbx_sifts_xref_db_num",
        "label_seq_id",
    ]

    table = pd.DataFrame(block.find("_atom_site.", loop_names), columns=loop_names)

    return table



def get_unp_residue_indices(model: pd.DataFrame, chain_id: str) -> set:

    # TODO: Consider using Gemmi's Selection class, although it might be slower
    # considering it uses the older MMDB parsing syntax.

    df = block_to_df(model)

    res_ids = df["pdbx_sifts_xref_db_num"][
        (df["label_asym_id"] == chain_id)
        & (df["group_PDB"] == "ATOM")
        & (df["pdbx_sifts_xref_db_num"] != "?")
    ].unique()

    # Sort and convert to int
    return set(map(int, res_ids))


def fill_missing_ends(
    all_data, all_res_ids, min_res_id, max_res_id, fill_type=np.full(3, np.nan)
):

    # Add NaNs to end of arrays
    for key, unp_res_id in all_res_ids.items():
        difference = max_res_id - max(unp_res_id)

        if difference == 0:
            continue

        appendage_array = np.full((difference, 3), np.nan)

        all_data[key] = np.concatenate((all_data[key], appendage_array))

    return all_data


def calculate_atom_distances(coordinates1, coordinates2):
    
    coord_difference = coordinates1 - coordinates2
    distances = np.linalg.norm(coord_difference, axis=1)

    return np.array(distances)







class PositionCalculator:
    def __init__(self, model: gemmi.Structure, unp_res_ids: set) -> None:
        self.model = model
        self.unp_res_ids = sorted(list(unp_res_ids))
        self.all_res_ids = range(1, max(unp_res_ids) + 1)

        # Calculate positions
        self.ca_coords = self.__extract_positions()

    def __extract_positions(self, atom_type="CA"):
        """
        :return: _description_
        :rtype: tuple( list(float), list(float) )
        """

        positions = []

        sorted_res_ids = sorted(list(self.unp_res_ids))

        index_position = 0  # Keep track of residue to extract
        for i in self.all_res_ids:
            if i in sorted_res_ids:
                ca_position = self.model[index_position][atom_type][0].pos
                ca_position = np.array([ca_position[0], ca_position[1], ca_position[2]])
                positions.append(ca_position)

                index_position += 1

            else:
                # Append NaNs
                positions += [np.array([np.nan, np.nan, np.nan])]

        return np.array(positions)

    def get_positions(self):
        return self.ca_coords
    




if __name__ == "__main__":


    # Configure logging 
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    # Collect arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--mmcif1", type=str, help="Path to mmCIF file")
    parser.add_argument("--mmcif2", type=str, help="Path to mmCIF file")
    parser.add_argument("--rt_matrices", type=str, help="Path to JSON file")
    parser.add_argument("--pdb_id1", type=str, help="PDB ID")
    parser.add_argument("--pdb_id2", type=str, help="PDB ID")
    parser.add_argument("--chain1", type=str, help="Chain ID")
    parser.add_argument("--chain2", type=str, help="Chain ID")
    args = parser.parse_args()

    # Entry labels
    entry1 = f"{args.pdb_id1}_{args.chain1}"
    entry2 = f"{args.pdb_id2}_{args.chain2}"

    # Load mmCIF files
    model1 = load_mmcif(args.mmcif1, return_model=True)
    model2 = load_mmcif(args.mmcif2, return_model=True)

    # Load RT matrices
    rt_matrices = load_RT_matrices(args.rt_matrices)
    rt_matrices1 = rt_matrices[f"{args.pdb_id1}_{args.chain1}"]["matrix"]
    rt_matrices2 = rt_matrices[f"{args.pdb_id2}_{args.chain2}"]["matrix"]

    # Convert RT matrices to Gemmi Transform objects
    rt_matrices1 = affine_to_transformation(rt_matrices1)
    rt_matrices2 = affine_to_transformation(rt_matrices2)
    
    # Apply RT matrices to models
    model1.transform_pos_and_adp(rt_matrices1)
    model2.transform_pos_and_adp(rt_matrices2)

    # Obtain residue indices
    block1 = load_mmcif(args.mmcif1, return_model=False)
    block2 = load_mmcif(args.mmcif2, return_model=False)

    unp_res_ids1 = get_unp_residue_indices(block1, args.chain1)
    unp_res_ids2 = get_unp_residue_indices(block2, args.chain2)
    all_res_ids = set()
    all_unp_res_ids = {}
    for i, j in zip((unp_res_ids1, unp_res_ids2), (entry1, entry2)):
        all_res_ids = all_res_ids.union(i)
        all_unp_res_ids[j] = i

    # Create atomic position objects
    all_atomic_coords = {}

    atomic_coords1 = PositionCalculator(model1[args.chain1], unp_res_ids1).get_positions()
    atomic_coords2 = PositionCalculator(model2[args.chain2], unp_res_ids2).get_positions()

    all_atomic_coords[entry1] = atomic_coords1
    all_atomic_coords[entry2] = atomic_coords2

    # Add missing residues to ends as NaN values
    all_coordinates = fill_missing_ends(
        all_atomic_coords,
        all_unp_res_ids,
        min(all_res_ids),
        max(all_res_ids),
        fill_type=np.array([np.nan, np.nan, np.nan]),
    )

    atomic_distances = calculate_atom_distances(
        all_coordinates[entry1], 
        all_coordinates[entry2]
    )

    fig, ax = plt.subplots(1,1, figsize=(10,5), tight_layout=True)

    ax.plot(range(0, max(all_res_ids)), atomic_distances, color="black", linewidth=0.5)
    ax.set_xlabel("UniProt Residue ID")
    ax.set_ylabel("Distance (\u212B)")
    ax.set_title(f"Euclidean distance between {entry1} vs {entry2}")

  
    # Save the plot to a file
    plot_filename = "atomic_distances_plot.png"
    plt.savefig(plot_filename)
    logging.info(f"Plot saved to {plot_filename}")

    # Convert atomic distances to a Pandas DataFrame
    df_data = {"UniProt Residue ID": range(0, max(all_res_ids)), "Distance (\u212B)": atomic_distances}
    df = pd.DataFrame(df_data)

    # Filter rows with Distance values between 0 and 2
    df_filtered = df[(df["Distance (\u212B)"] >= 0) & (df["Distance (\u212B)"] <= 2)]

    # Save the DataFrame to an Excel file
    excel_filename = "atomic_distances_table_filtered.xlsx"
    df_filtered.to_excel(excel_filename, index=False)
    logging.info(f"Filtered Atomic distances table (0-2 \u212B) saved to {excel_filename}")

    # Save the DataFrame to a CSV file
    csv_filename = "atomic_distances_table_filtered.csv"
    df_filtered.to_csv(csv_filename, index=False)
    logging.info(f"Filtered Atomic distances table (0-2 \u212B) saved to {csv_filename}")


    # Save the DataFrame to an Excel file

    excel_filename = "atomic_distances_table.xlsx"
    try:
        df.to_excel(excel_filename, index=False)
        logging.info(f"Atomic distances table saved to {excel_filename}")
    except Exception as e:
        logging.error(f"Error saving Excel file: {e}")

    # Save the DataFrame to a CSV file (added line)
    csv_filename = "atomic_distances_table.csv"
    df.to_csv(csv_filename, index=False)
    logging.info(f"Atomic distances table saved to {csv_filename}")
