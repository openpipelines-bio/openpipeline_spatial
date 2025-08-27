import fnmatch
import zipfile
import tempfile
from pathlib import Path
from typing import Union


def unzip_archived_folder(
  archived_folder: Union[str, Path]
) -> Union[str, Path]:
    """
    Extracts a ZIP archive to a temporary directory and returns the path to the extracted folder.

    Args:
        zip_path (Union[str, Path]): Path to the ZIP archive.

    Returns:
        extracted_path (Union[str, Path]): Path to the extracted folder inside the temporary directory.
    """

    temp_dir = Path(tempfile.TemporaryDirectory().name)
    with zipfile.ZipFile(archived_folder, "r") as archive:
        archive.extractall(temp_dir)

    return temp_dir


def extract_selected_files_from_zip(
    zip_path: Union[str, Path],
    members: list[Union[str, Path]]
) -> Union[str, Path]:
    """
    Extracts selected files (supports glob patterns) from a ZIP archive to a temporary directory.

    Args:
        zip_path (Union[str, Path]): Path to the ZIP archive.
        members (list[str]): List of file paths within the archive to extract.

    Returns:
        Path: Path to the extraction directory.
    """

    temp_dir = Path(tempfile.TemporaryDirectory().name)

    with zipfile.ZipFile(zip_path, "r") as archive:
        all_files = archive.namelist()
        selected = set()
        for pattern in members:
            selected.update(fnmatch.filter(all_files, str(pattern)))
        for member in selected:
            archive.extract(member, temp_dir)

    return temp_dir
