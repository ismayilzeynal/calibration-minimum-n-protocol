"""Copy (read-only) the small Eedi metadata members needed by the revision
analyses into revision-2026-11/data/. Usage:
    python extract_metadata.py <path to Eedi data.zip>
"""
import hashlib
import json
import os
import sys
import zipfile

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
MEMBERS = ["question_metadata_task_1_2.csv", "question_metadata_task_3_4.csv",
           "student_metadata_task_1_2.csv", "student_metadata_task_3_4.csv",
           "subject_metadata.csv"]


def main(zip_path):
    os.makedirs(OUT, exist_ok=True)
    z = zipfile.ZipFile(zip_path)
    log = {"source_zip": os.path.basename(zip_path),
           "source_zip_bytes": os.path.getsize(zip_path), "members": {}}
    for m in MEMBERS:
        b = z.read("data/metadata/" + m)
        with open(os.path.join(OUT, m), "wb") as fh:
            fh.write(b)
        log["members"][m] = {"bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()}
        print(m, len(b), log["members"][m]["sha256"])
    with open(os.path.join(OUT, "METADATA_SOURCE.json"), "w") as fh:
        json.dump(log, fh, indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
