import pandas as pd

def main():
    src_data_path = "/Volumes/Elements/vastai/ppmi/ppmi_20260415_170143/data.parquet"

    data = pd.read_parquet(src_data_path)
    print(data.head())
    data["Sample_Group"].value_counts()
    
if __name__ == "__main__":
    main()