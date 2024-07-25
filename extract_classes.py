from Bio import SeqIO
import pandas as pd

# Load the GenBank file
gbk_file = "umap/GCA_020531285_antiSMASH/JAJBGS010000009.region001.gbk"

# Initialize lists to store extracted information
file_names = []
protocluster_counts = []
protocluster_classes = []

# Parse the GenBank file
for record in SeqIO.parse(gbk_file, "genbank"):
    file_names.append(record.id)
    
    # Extract protocluster information
    protoclusters = record.features_of_type("protocluster")
    protocluster_counts.append(len(protoclusters))
    
    # Extract protocluster classes
    protocluster_class = [protocluster.qualifiers.get("category")[0] for protocluster in protoclusters]
    protocluster_classes.append(protocluster_class)

# Create a dataframe
data = {
    "File Name": file_names,
    "Protocluster Count": protocluster_counts,
    "Protocluster Classes": protocluster_classes
}

df = pd.DataFrame(data)

# Display the dataframe
print(df)
