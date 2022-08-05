from ihkapy.fileio.metadata_io import parse_metadata
import os 

if __name__ == "__main__":
    rel_data_path = "./test_data"
    session_basename = "AC75a-5 DOB 072519_TS_2020-03-23_17_30_04"
    session_metadata_path = os.path.join(rel_data_path,session_basename + ".txt")
    fm,seiz_info = parse_metadata(session_metadata_path)

    print("Frontmatter")
    print(fm)

    print("\n\nSeizure info")
    print(seiz_info)


