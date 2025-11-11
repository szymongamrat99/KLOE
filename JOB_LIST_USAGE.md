# Job List File Usage

## Overview

The program now supports reading file lists from a text file passed as a command-line argument. This allows you to specify exactly which ROOT files to analyze without relying on the standard directory scanning method.

## File Format

### Filename Format
The job list file must follow this naming convention:
```
job_v{version}_{type}_{luminosity}_inv_pb_{number}.txt
```

**Parameters:**
- `{version}`: Analysis version number (e.g., 1, 2, 3)
- `{type}`: Analysis type - one of: `data`, `all_phys`, `all_phys2`, `all_phys3`
- `{luminosity}`: Luminosity value in pb^-1 (e.g., 5000, 10000)
- `{number}`: Job number (e.g., 001, 002)

**Examples of valid filenames:**
```
job_v1_data_5000_inv_pb_001.txt
job_v2_all_phys_10000_inv_pb_001.txt
job_v1_all_phys3_2500_inv_pb_042.txt
```

### File Content
The file should contain absolute paths to ROOT files, one per line:

```
/storage/gpfs_data/kloe/dataroot/DBV-26/DK0/prod2root_dk0_30548.root
/storage/gpfs_data/kloe/dataroot/DBV-26/DK0/prod2root_dk0_30549.root
/storage/gpfs_data/kloe/dataroot/DBV-26/DK0/prod2root_dk0_30550.root
/storage/gpfs_data/kloe/dataroot/DBV-26/DK0/prod2root_dk0_30551.root
```

**Notes:**
- Empty lines are ignored
- Lines starting with `#` are treated as comments and ignored
- Whitespace at the beginning and end of lines is trimmed
- Each path must point to a valid `.root` file

## Usage

### Command Line

Run the program with the job list file as an argument:

```bash
./KLSPM00 job_v1_data_5000_inv_pb_001.txt
```

### What Happens

1. The program validates the filename format
2. Reads all file paths from the file
3. Verifies that each file exists and has `.root` extension
4. Adds all files to the TChain
5. Automatically determines the analysis type from the filename (data, all_phys, etc.)
6. Performs the initial analysis using the specified files

### Example Workflow

```bash
# Create a job list file
cat > job_v1_data_5000_inv_pb_001.txt << 'EOF'
/path/to/file1.root
/path/to/file2.root
/path/to/file3.root
EOF

# Run the analysis
./KLSPM00 job_v1_data_5000_inv_pb_001.txt
```

## Error Handling

The program will report errors in the following cases:

1. **File not found**: The job list file doesn't exist
   ```
   ERROR: Job list file does not exist: job_v1_data_5000_inv_pb_001.txt
   ```

2. **Invalid filename format**: The filename doesn't match the expected pattern
   ```
   ERROR: Job list file has invalid name format.
   Expected format: job_v{version}_{type}_{luminosity}_inv_pb_{number}.txt
   Valid types: data, all_phys, all_phys2, all_phys3
   ```

3. **Missing files in list**: A path in the file doesn't exist
   ```
   WARNING: File does not exist: /path/to/missing_file.root
   ```

4. **Non-ROOT files**: A file in the list doesn't have `.root` extension
   ```
   WARNING: File is not ROOT file: /path/to/file.txt
   ```

5. **Empty file list**: The file contains no valid file paths
   ```
   ERROR: No files found in job list file: job_v1_data_5000_inv_pb_001.txt
   ```

## Logging

All operations are logged in the program's log file. The log will contain:
- Number of files loaded from the job list
- Analysis type determined from filename
- List of all files added to the TChain
- Total number of entries loaded

## Standard Operation (Without Job List)

If no argument is provided, the program operates in standard mode:

```bash
./KLSPM00
```

In this case:
- The program uses `initialAnalysisExecution` flag from config
- Or presents interactive menu for file selection
- Standard directory scanning is used to find files
