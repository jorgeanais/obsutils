from pathlib import Path

from astropy.coordinates import SkyCoord
import astropy.units as u


def pointing_error(x, y):
    """Compute the difference in arcseconds between the input coordinates and the center of the image for Swope C3 quadrant"""
    X_CENTER = 1024
    Y_CENTER = 1028
    SCALE = 0.435
    print(f"x offset:  {(x - X_CENTER)*SCALE:.1f}")
    print(f"x offset:  {(y - Y_CENTER)*SCALE:.1f}")


def apply_offset(ra, dec, ra_offset=375, dec_offset=425):
    """
    Apply offsets to Right Ascension (RA) and Declination (DEC) coordinates.

    Parameters:
    ra (str or float): Object right ascension.
    dec (str or float): Object declination.
    ra_offset (float): Offset in arcseconds to be added to the RA coordinate (default is 375).
    dec_offset (float): Offset in arcseconds to be added to the DEC coordinate (default is 425).

    Returns:
    tuple: Tuple containing the modified RA and DEC coordinates.
    """
    # Convert coordinates to SkyCoord object
    initial_coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))

    # Apply offsets
    new_ra = initial_coord.ra + ra_offset * u.arcsec
    new_dec = initial_coord.dec + dec_offset * u.arcsec

    print(
        f"RA = {initial_coord.ra.to_string(unit=u.hourangle, decimal=True)} {new_ra.to_string(unit=u.hourangle, decimal=True)}"
    )  # TODO: Check if hourangle or deg
    print(
        f"DEC = {initial_coord.dec.to_string(unit=u.deg, decimal=True)} {new_dec.to_string(unit=u.deg, decimal=True)}"
    )
    return new_ra, new_dec


def parse_input_file(path: Path) -> list[str]:
    """
    Parse the input file containing the coordinates of the objects to be observed.

    Parameters:
    path (Path): Path object pointing to the input file.

    Returns:
    list: List containing the name and coordinates of the objects to be observed with the offset applied.
    """

    with open(path, "r") as f:
        lines = f.readlines()
        non_empty_lines = [line.strip() for line in lines if line.strip()]

    object_list = [line.split("\t") for line in non_empty_lines]

    output = []
    for object_ in object_list:
        shifted_ra, shifted_dec = apply_offset(object_[1], object_[2])
        output.append(
            f"%{object_[0]}%{shifted_ra.to_string(unit=u.hourangle, decimal=True)}%{shifted_dec.to_string(unit=u.deg)}"
        )

    return output


def save_to_output_file(path: Path, content_list: list[str]) -> None:
    """
    Save the results of parse_input_file to an output file

    Parameters:
    path (Path): Path object pointing to the input file.
    content_list (list): List containing the name and coordinates of the objects as from parse_input_file.
    """
    with open(path, "w") as f:
        for item in content_list:
            # Write each item to the file
            f.write(item + "\n")


if __name__ == "__main__":
    # print(pointing_error(1189, 1045))
    # new_ra, new_dec = apply_offset("7 52 42.6", "14 50 21.4")
    # print(new_ra, new_dec)
    l = parse_input_file("data/input.txt")
    save_to_output_file("data/output.txt", l)
