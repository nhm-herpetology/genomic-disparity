This protocol was developed as part of the ULTRAMOD project https://github.com/AshwiniVM/ULTRAMOD

# Step 1: Curating chromosome-level genome assemblies
We will need to download chromosome-level genome assemblies from NCBI or other repository. 

# Step 2: Selecting landmark sequences and formatting for reference mapping
Landmarks can be any conserved sequence that can be aligned to genomes in your dataset, but we have used ultraconserved elements (UCEs) as an example.  

# Step 3: Reference mapping of landmarks to genomes
BWA is used to map the landmarks to different chromosomes

# Step 4: Identify those chromosomes that contain very similar landmarks

# Step 5: Remove non-common landmarks and curate chromosome sets
