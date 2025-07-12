
<!-- README.md is generated from README.Rmd. Please edit that file -->

# 🧪 labtools

<div align="center">

<img src="https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white" alt="R" height="32"/>

## ⚗️ Analytical Chemistry Laboratory Data Processing Toolkit

*A comprehensive R package for streamlining analytical chemistry
workflows*

<!-- Badges -->

<img src="https://img.shields.io/github/workflow/status/QizhiSu/labtools/R-CMD-check?style=flat-square&logo=github&logoColor=white&labelColor=21262d&color=3fb950" alt="R-CMD-check"/>
<img src="https://img.shields.io/badge/License-MIT-7c3aed?style=flat-square&labelColor=21262d" alt="License"/>
<img src="https://img.shields.io/badge/CRAN-not%20yet-ff7b72?style=flat-square&labelColor=21262d" alt="CRAN"/>
<img src="https://img.shields.io/github/stars/QizhiSu/labtools?style=flat-square&logo=github&logoColor=white&labelColor=21262d&color=d29922" alt="Stars"/>
<img src="https://img.shields.io/github/issues/QizhiSu/labtools?style=flat-square&logo=github&logoColor=white&labelColor=21262d&color=58a6ff" alt="Issues"/>
<img src="https://img.shields.io/github/last-commit/QizhiSu/labtools?style=flat-square&logo=github&logoColor=white&labelColor=21262d&color=7c3aed" alt="Last Commit"/>

<!-- Language Toggle -->

[🇺🇸 English](#english) \| [🇨🇳 中文](#中文)

</div>

------------------------------------------------------------------------

## 🎯 **About This Package**

> **labtools** is specifically designed for our lab’s **analytical
> chemistry workflows** and daily data processing needs. This package
> prioritizes **ease of use for researchers who may not be R experts**
> over maximum flexibility.
>
> **Design Philosophy:** - **Laboratory-focused**: Built around common
> analytical chemistry tasks and workflows - **User-friendly**: Simple
> function calls with sensible defaults for non-R experts -
> **Workflow-oriented**: Functions designed to work together in typical
> lab data processing pipelines - **Trade-offs**: Some flexibility is
> sacrificed for simplicity and ease of use
>
> **Please use according to your needs** - if you require maximum
> flexibility, consider using the underlying packages directly (webchem,
> rcdk, etc.).

------------------------------------------------------------------------

## <a id="english"></a>🌟 **Key Features**

| Feature                              | Description                                                                          |
|--------------------------------------|--------------------------------------------------------------------------------------|
| 🔍 **Chemical Metadata Extraction**  | Enhanced automated retrieval from PubChem with robust error handling and CAS parsing |
| 🧬 **Structure Database Management** | Build and visualize chemical structure databases with advanced filtering             |
| 📊 **MS Data Processing**            | Export optimized databases for MS-FINDER and MS-DIAL software                        |
| 🔬 **2D GC-MS Analysis**             | Process Canvas exports and combine multi-sample data with precision                  |
| ⚗️ **Chemical Structure Conversion** | SMILES to MOL file conversion with 2D coordinate generation                          |
| 📈 **Semi-quantification Tools**     | Advanced analytical quantification with machine learning integration                 |
| 🎯 **Interactive Visualization**     | Shiny-based chemical structure navigation and spectrum plotting                      |

> **💡 New to R?** For detailed help on any function, use
> `?function_name` in R console (e.g., `?extract_cid`)

### 🆕 Latest Updates (v0.3.00)

- **Enhanced metadata extraction**: Improved `parse_cas_clean()`,
  `extract_cid()`, and `extract_meta()` functions
- **Better dependency management**: `mspcompiler` is now optional -
  install only when needed
- **Robust error handling**: Comprehensive input validation and
  informative error messages
- **Improved documentation**: Detailed examples and parameter
  explanations for all functions
- **Enhanced synonyms extraction**: Fixed and optimized PubChem synonyms
  retrieval
- **ASCII compliance**: All code now uses ASCII characters for better
  portability
- **Comprehensive testing**: All functions pass R CMD check with zero
  errors and warnings

------------------------------------------------------------------------

## 🚀 **Installation**

<div class="code-block">

### Development Version (Recommended)

``` r
# Install devtools if you haven't already
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install labtools from GitHub
devtools::install_github("QizhiSu/labtools")
```

### System Requirements

- **R** ≥ 4.0.0
- **Java** ≥ 8 (for rcdk package)
- **Internet connection** (for PubChem API access)

### Optional Dependencies

Some functions require additional packages that are not automatically
installed:

- **mspcompiler**: Required for MSP library filtering (`filter_msp()`
  function)

  ``` r
  # Install when needed
  devtools::install_github("QizhiSu/mspcompiler")
  ```

### Getting Help

- **Function help**: Use `?function_name` (e.g., `?extract_cid`)
- **Package help**: Use `help(package = "labtools")`
- **Examples**: Use `example(function_name)` to run examples

</div>

------------------------------------------------------------------------

## ⚡ **Quick Start**

<div class="code-block">

``` r
library(labtools)
library(dplyr)

# 🔍 1. Extract chemical metadata from PubChem
data <- data.frame(
  Name = c("Caffeine", "Aspirin", "Glucose"),
  CAS = c("58-08-2", "50-78-2", "50-99-7")
)

# Extract CIDs and comprehensive metadata
data_with_cid <- extract_cid(data, name_col = "Name", cas_col = "CAS")
data_complete <- extract_meta(data_with_cid, cas = TRUE, flavornet = TRUE)

# 📊 2. Export for MS software
export4msdial(data_complete, polarity = "pos", output = "database_msdial.txt")
export4msfinder(data_complete, output = "database_msfinder.txt")

# ⚗️ 3. Convert SMILES to MOL files
smiles_data <- data.frame(
  ID = c("Caffeine", "Aspirin"),
  SMILES = c("Cn1cnc2c1c(=O)n(c(=O)n2C)C", "CC(=O)OC1=CC=CC=C1C(=O)O")
)
export_smiles_to_mol(smiles_data, output_dir = "mol_files")

# 🎯 4. Interactive chemical structure browser
navigate_chem(data_complete)  # Opens Shiny app
```

<div class="info-box">

<strong>💡 Tip:</strong> Use <code>?extract_cid</code> or
<code>?export_smiles_to_mol</code> to see detailed parameter
descriptions and more examples!

</div>

</div>

------------------------------------------------------------------------

## 📚 **Detailed Function Examples**

### 🔍 Chemical Metadata Extraction

<div class="code-block">

``` r
library(labtools)

# Extract CIDs from chemical identifiers
compounds <- data.frame(
  Name = c("Caffeine", "Theobromine"),
  CAS = c("58-08-2", "83-67-0")
)

# Extract CIDs with proper error handling
# Use ?extract_cid for full parameter details
compounds_cid <- extract_cid(
  data = compounds,           # Input data frame
  name_col = "Name",         # Column with chemical names (default: "Name")
  cas_col = "CAS",           # Column with CAS numbers (default: "CAS")
  inchikey_col = "InChIKey", # Column with InChIKeys (default: "InChIKey")
  verbose = TRUE             # Show progress messages (default: TRUE)
)

# Extract comprehensive metadata
# Use ?extract_meta for all options
compounds_meta <- extract_meta(
  data = compounds_cid,      # Data frame with CID column (required)
  cas = TRUE,                # Extract CAS numbers (default: FALSE)
  flavornet = TRUE,          # Extract Flavornet data (default: FALSE)
  synonyms = TRUE,           # Extract synonyms (default: FALSE)
  uses = TRUE,               # Extract compound uses (default: FALSE)
  verbose = TRUE             # Show progress (default: TRUE)
)

# Add chemical classification using ClassyFire
# Use ?extract_classyfire for details
compounds_classified <- extract_classyfire(
  data = compounds_meta,     # Data frame with InChIKey column
  inchikey_col = "InChIKey", # InChIKey column name (default)
  name_col = "Name"          # Name column for progress (default)
)
```

**Key Parameters Explained:** - `name_col`, `cas_col`, `inchikey_col`:
Specify which columns contain identifiers (defaults: “Name”, “CAS”,
“InChIKey”) - `verbose`: Set to `FALSE` to suppress progress messages
(default: TRUE) - `cas`, `flavornet`, `synonyms`, `uses`: Enable
extraction of specific metadata types (all default: FALSE) - Functions
automatically handle network errors and provide informative progress
messages

</div>

### 📤 Database Export

<div class="code-block">

``` r
# Export for MS-FINDER software
# Use ?export4msfinder for details
export4msfinder(
  data = compounds_classified,              # Data with required columns
  output = "msfinder_database.txt"          # Output file path (default)
)

# Export for MS-DIAL (positive mode)
# Use ?export4msdial for all options
export4msdial(
  data = compounds_classified,              # Data with required columns
  polarity = "pos",                         # ESI polarity: "pos" or "neg" (default: "pos")
  output = "msdial_pos.txt"                # Output file path (default)
)

# Export for MS-DIAL (negative mode)
export4msdial(
  data = compounds_classified,
  polarity = "neg",                         # Negative mode for different adducts
  output = "msdial_neg.txt"
)
```

**Required Columns:** - **MS-FINDER**: Name, InChIKey, CID, ExactMass,
Formula, SMILES - **MS-DIAL**: Name, ExactMass, SMILES, InChIKey

**Polarity Options:** - `"pos"`: Generates \[M+H\]+, \[M+Na\]+, \[M+K\]+
adducts - `"neg"`: Generates \[M-H\]-, \[M+Cl\]-, \[M+HCOO\]- adducts

</div>

### ⚗️ SMILES to MOL Conversion

<div class="code-block">

``` r
# Prepare SMILES data
smiles_data <- data.frame(
  Compound_ID = c("CAFF_001", "THEO_002"),
  SMILES_String = c(
    "Cn1cnc2c1c(=O)n(c(=O)n2C)C",      # Caffeine
    "Cn1cnc2c1c(=O)[nH]c(=O)n2C"       # Theobromine
  ),
  stringsAsFactors = FALSE
)

# Convert to MOL files with 2D coordinates
# Use ?export_smiles_to_mol for all parameters
result <- export_smiles_to_mol(
  df = smiles_data,                     # Data frame with ID and SMILES columns
  id_col = "Compound_ID",               # Column with compound IDs (default: "ID")
  smiles_col = "SMILES_String",         # Column with SMILES strings (default: "SMILES")
  output_dir = "molecular_structures"   # Output directory (default: "mol_files")
)

# Check conversion results
print(result)
# Returns: list(success = 2, failed = 0, skipped = 0)
```

**Function Features:** - **Input validation**: Checks for valid data
frame and required columns - **Error handling**: Continues processing
even if some SMILES fail - **2D coordinates**: Automatically generates
2D coordinates for visualization - **Progress tracking**: Returns
summary of successful/failed conversions - **Robust processing**:
Handles NA/NULL values gracefully

</div>

### 🔬 Canvas Data Processing

<div class="code-block">

``` r
# Process Canvas 2D GC-MS data with detailed parameters
# Use ?read_canvas for complete parameter documentation
canvas_data <- read_canvas(
  path = "path/to/canvas/files",        # Directory containing Canvas .txt files
  ri_iden_tol = 30,                     # RI tolerance for identification (default: 30)
  ri_p_iden_tol = 100,                  # Predicted RI tolerance for identification (default: 100)
  ri_align_tol = 50,                    # RI tolerance for alignment (default: 50)
  rt_2d_align_tol = 0.1,                # 2D RT tolerance for alignment (default: 0.1)
  keep = "area"                         # Data to keep: "area", "height", or "both" (default: "area")
)

# Normalize by internal standard (D8)
# Use ?normalize_area for details
normalized_data <- normalize_area(
  df = canvas_data,                     # Data frame with peak area data
  start_col = 12                        # Starting column index of area data (default: 12)
)

# Filter areas based on sample codes
# Use ?keep_area for filtering options
filtered_data <- keep_area(
  df = normalized_data,                 # Data frame with peak area data
  sam_code = sample_codes,              # Data frame mapping sample codes to names
  start_col = 12,                       # Starting column index (default: 12)
  keep_bk = TRUE,                       # Keep blank samples (default: TRUE)
  keep_d8 = TRUE                        # Keep D8 internal standard (default: TRUE)
)
```

**Key Parameters Explained:** - `ri_iden_tol`, `ri_p_iden_tol`:
Retention index tolerances for compound identification - `ri_align_tol`,
`rt_2d_align_tol`: Tolerances for aligning peaks across samples -
`keep`: Determines which peak data to retain (area is most common for
quantification) - `start_col`: Column index where peak area data begins
(typically after metadata columns) - `keep_bk`, `keep_d8`: Control
whether to retain blank samples and internal standards

</div>

### 📈 Semi-quantification Analysis

<div class="code-block">

``` r
# Select appropriate standards using molecular similarity
# Use ?select_std for detailed parameter explanation
standards_assigned <- select_std(
  std_md = reference_standards,         # Standards data frame with response and descriptors
  std_res_col = 3,                      # Column index of response variable (or FALSE for all descriptors)
  std_md1_col = 14,                     # Starting column index of molecular descriptors in standards
  data_md = target_compounds,           # Target compounds data frame with descriptors
  data_md1_col = 5,                     # Starting column index of descriptors in targets
  top_npct_md = 30                      # Percentage of top molecular descriptors to use (default: 20)
)

# Calculate concentrations using assigned standards
# Use ?calculate_con for concentration calculation details
concentrations <- calculate_con(
  df = standards_assigned,              # Data frame after select_std() processing
  sam_weight = sample_weights,          # Data frame with sample weights and volumes (can be NULL)
  start_col = 12                        # Starting column index of peak area data (default: 12)
)

# Organize and summarize concentration data
# Use ?organize_con for result organization options
final_results <- organize_con(
  df = concentrations,                  # Data frame after calculate_con() processing
  sam_code = sample_codes,              # Data frame with sample codes and names
  start_col = 12,                       # Starting column index of concentration data (default: 12)
  digits = 3,                           # Decimal places for rounding (default: 2)
  na2zero = TRUE,                       # Replace NA with zero for calculations (default: TRUE)
  bind_mean_sd = TRUE                   # Combine mean and SD in single columns (default: TRUE)
)
```

**Semi-quantification Workflow:** 1. **Standard Selection**: Uses
molecular similarity to assign appropriate calibration standards 2.
**Concentration Calculation**: Applies standard curves to calculate
concentrations 3. **Result Organization**: Summarizes data with
statistics and proper formatting

**Key Parameters:** - `std_res_col`: Set to `FALSE` to use all molecular
descriptors without feature selection - `top_npct_md`: Controls how many
molecular descriptors to use (higher = more descriptors) - `sam_weight`:
Include sample weights for concentration per gram calculations -
`bind_mean_sd`: Combines mean±SD into single readable columns

</div>

### 📊 Spectral Analysis

<div class="code-block">

``` r
# Convert spectrum string to numeric vector
# Use ?update_spectrum for format details
spectrum_str <- "100:30 120:80 145:60 170:90 200:20"
spectrum <- update_spectrum(
  spectrum_str = spectrum_str,          # Spectrum in "mz:intensity mz:intensity" format
  start_mz = 50,                        # Starting m/z value (default: 50)
  end_mz = 500,                         # Ending m/z value (default: 500)
  mz_step = 1,                          # m/z step size (default: 1)
  digits = 0                            # Decimal places for m/z rounding (default: 0)
)

# Create interactive spectrum plot
# Use ?plot_spectrum for plotting options
plot_spectrum(
  spectrum = spectrum,                  # Named numeric vector (m/z as names, intensity as values)
  range = 10,                           # Bin width for peak labeling (default: 10)
  threshold = 5,                        # Minimum intensity (%) for peak labeling (default: 1)
  max_ticks = 20                        # Maximum number of x-axis ticks (default: 20)
)

# Calculate spectral similarity using cosine similarity
# Use ?cosine_similarity for similarity calculation details
spec1 <- update_spectrum("100:30 120:80 145:60")
spec2 <- update_spectrum("100:25 120:85 145:55")
similarity <- cosine_similarity(
  spectrum1 = spec1,                    # First spectrum (named numeric vector)
  spectrum2 = spec2                     # Second spectrum (named numeric vector)
)

# Convert between MSP and data frame formats
# Use ?msp2df and ?df2msp for format conversion
msp_data <- msp2df(msp_file_path)      # Convert MSP file to data frame
df2msp(df_data, "output.msp")          # Convert data frame to MSP file
```

**Spectral Data Formats:** - **Input**: Spectrum strings in
“mz:intensity mz:intensity” format - **Processing**: Converts to named
numeric vectors for analysis - **Output**: Interactive plots and
similarity scores

**Key Functions:** - `update_spectrum()`: Converts text format to
numeric vectors - `plot_spectrum()`: Creates interactive mass spectrum
plots - `cosine_similarity()`: Calculates spectral similarity (0-1
scale) - `msp2df()` / `df2msp()`: Convert between MSP files and data
frames

</div>

### 🔍 MSP Library Filtering

<div class="code-block">

``` r
# Note: MSP filtering requires the mspcompiler package
# Install if needed: devtools::install_github("QizhiSu/mspcompiler")

# Prepare compound list for filtering
compounds_to_keep <- data.frame(
  Name = c("Caffeine", "Aspirin", "Glucose"),
  InChIKey = c("RYYVLZVUVIJVGH-UHFFFAOYSA-N",
               "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
               "WQZGKKKJIJFFOK-GASJEMHNSA-N"),
  stringsAsFactors = FALSE
)

# Filter MSP spectral library based on compound list
# Use ?filter_msp for filtering options
filter_msp(
  msp = "NIST_library.msp",             # Path to the MSP library file
  cmp_list = compounds_to_keep,         # Data frame with compounds to keep (Name and InChIKey columns required)
  keep_napd8 = TRUE,                    # Add Naphthalene-D8 to the filter list (default: TRUE)
  output = "filtered_library.msp"       # Output path for filtered MSP file
)

# Read MS-DIAL output files with filtering options
# Use ?read_msdial for reading options
msdial_data <- read_msdial(
  file_path = "msdial_results.txt",     # Path to MS-DIAL output file
  keep_unknown = FALSE,                 # Keep unknown compounds (default: FALSE)
  keep_spectrum = TRUE,                 # Keep spectrum data (default: FALSE)
  keep_mean_sd = FALSE                  # Keep mean and SD columns (default: FALSE)
)
```

**MSP Library Management:** - **Purpose**: Filter large spectral
libraries to contain only compounds of interest - **Requirements**:
Compound list must have Name and InChIKey columns - **Features**:
Automatically includes Naphthalene-D8 internal standard - **Output**:
Filtered MSP file ready for MS software

**MS-DIAL Integration:** - **Input**: MS-DIAL alignment result files -
**Options**: Control retention of unknown compounds and spectral data -
**Output**: Clean data frame ready for further analysis

</div>

------------------------------------------------------------------------

## 🔬 **Advanced Workflows**

### Complete Database Creation Workflow

``` r
# Step 1: Prepare data
chemicals <- data.frame(
  Name = c("Benzene", "Toluene", "Xylene"),
  CAS = c("71-43-2", "108-88-3", "1330-20-7")
)

# Step 2: Extract comprehensive metadata
chemicals_cid <- extract_cid(chemicals, name_col = "Name", cas_col = "CAS")
chemicals_meta <- extract_meta(chemicals_cid, cas = TRUE, synonyms = TRUE)
chemicals_class <- extract_classyfire(chemicals_meta)

# Step 3: Export for different platforms
export4msdial(chemicals_class, polarity = "pos", output = "pos_database.txt")
export4msdial(chemicals_class, polarity = "neg", output = "neg_database.txt")
export4msfinder(chemicals_class, output = "msfinder_database.txt")

# Step 4: Generate structure files
export_smiles_to_mol(chemicals_class, output_dir = "structures")
```

### 2D GC-MS Processing Pipeline

``` r
# Process Canvas data
canvas_data <- read_canvas("canvas_files", keep = "area")

# Extract metadata
canvas_meta <- extract_cid(canvas_data, name_col = "Name", cas_col = "CAS")
canvas_complete <- extract_meta(canvas_meta, cas = TRUE)

# Normalize and analyze
canvas_normalized <- normalize_area(canvas_complete)
canvas_filtered <- keep_area(canvas_normalized, sample_codes)

# Semi-quantification
standards <- select_std(ref_standards, 3, 14, canvas_filtered, 5)
concentrations <- calculate_con(standards, sample_weights)
results <- organize_con(concentrations, sample_codes)
```

------------------------------------------------------------------------

## <a id="中文"></a>🧪 **labtools**: 分析化学实验室数据处理工具包

> **专为分析化学工作流程设计的综合性R包**

<div class="info-box">

**labtools**
专门为**分析化学实验室工作流程**和日常数据处理需求而设计。本包优先考虑**不太熟悉R语言的研究人员的易用性**，而非最大的灵活性。

**设计理念：** - **实验室导向**：围绕常见的分析化学任务和工作流程构建 -
**用户友好**：为非R专家提供简单的函数调用和合理的默认设置 -
**工作流导向**：函数设计为在典型的实验室数据处理流程中协同工作 -
**权衡取舍**：为了简单性和易用性，牺牲了一些灵活性

**请根据您的需求使用** -
如果您需要最大的灵活性，请考虑直接使用底层包（webchem、rcdk等）。

</div>

<div class="feature-table">

### 🌟 **主要功能**

| 功能                         | 描述                                      |
|------------------------------|-------------------------------------------|
| 🔍 **化学元数据提取**        | 从PubChem数据库自动检索，具有智能重试机制 |
| 🧬 **结构数据库管理**        | 构建和可视化化学结构数据库，支持高级过滤  |
| 📊 **质谱数据处理**          | 为MS-FINDER和MS-DIAL软件导出优化数据库    |
| 🔬 **二维气相色谱-质谱分析** | 精确处理Canvas导出数据并合并多样本数据    |
| ⚗️ **化学结构转换**          | SMILES到MOL文件转换，支持2D坐标生成       |
| 📈 **半定量工具**            | 集成机器学习的高级分析定量方法            |
| 🎯 **交互式可视化**          | 基于Shiny的化学结构导航和光谱绘图         |

</div>

<div class="info-box">

<strong>💡 R语言新手？</strong>
对于任何函数的详细帮助，请在R控制台中使用 <code>?函数名</code>
（例如：<code>?extract_cid</code>）

</div>

### 🚀 **安装**

<div class="code-block">

``` r
# 安装devtools（如果尚未安装）
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# 从GitHub安装labtools
devtools::install_github("QizhiSu/labtools")
```

### 系统要求

- **R** ≥ 4.0.0
- **Java** ≥ 8 (rcdk包需要)
- **网络连接** (访问PubChem API)

### 获取帮助

- **函数帮助**: 使用 `?函数名` (例如：`?extract_cid`)
- **包帮助**: 使用 `help(package = "labtools")`
- **示例**: 使用 `example(函数名)` 运行示例

</div>

### ⚡ **快速开始**

<div class="code-block">

``` r
library(labtools)
library(dplyr)

# 🔍 1. 从PubChem提取化学元数据
data <- data.frame(
  Name = c("咖啡因", "阿司匹林", "葡萄糖"),
  CAS = c("58-08-2", "50-78-2", "50-99-7")
)

# 提取CID和综合元数据
data_with_cid <- extract_cid(data, name_col = "Name", cas_col = "CAS")
data_complete <- extract_meta(data_with_cid, cas = TRUE, flavornet = TRUE)

# 📊 2. 导出用于质谱软件
export4msdial(data_complete, polarity = "pos", output = "database_msdial.txt")
export4msfinder(data_complete, output = "database_msfinder.txt")

# ⚗️ 3. 转换SMILES为MOL文件
smiles_data <- data.frame(
  ID = c("咖啡因", "阿司匹林"),
  SMILES = c("Cn1cnc2c1c(=O)n(c(=O)n2C)C", "CC(=O)OC1=CC=CC=C1C(=O)O")
)
export_smiles_to_mol(smiles_data, output_dir = "mol_files")

# 🎯 4. 交互式化学结构浏览器
navigate_chem(data_complete)  # 打开Shiny应用
```

<div class="info-box">

<strong>💡 提示:</strong> 使用 <code>?extract_cid</code> 或
<code>?export_smiles_to_mol</code> 查看详细的参数说明和更多示例！

</div>

</div>

### 📚 **详细函数示例**

#### 🔍 化学元数据提取

<div class="code-block">

``` r
library(labtools)

# 从化学标识符提取CID
compounds <- data.frame(
  Name = c("咖啡因", "可可碱"),
  CAS = c("58-08-2", "83-67-0")
)

# 提取CID并进行适当的错误处理
# 使用 ?extract_cid 查看完整参数详情
compounds_cid <- extract_cid(
  data = compounds,           # 输入数据框
  name_col = "Name",         # 包含化学名称的列 (默认: "Name")
  cas_col = "CAS",           # 包含CAS号的列 (默认: "CAS")
  inchikey_col = "InChIKey", # 包含InChIKey的列 (默认: "InChIKey")
  verbose = TRUE             # 显示进度信息 (默认: TRUE)
)

# 提取综合元数据
# 使用 ?extract_meta 查看所有选项
compounds_meta <- extract_meta(
  data = compounds_cid,      # 包含CID列的数据框 (必需)
  cas = TRUE,                # 提取CAS号 (默认: FALSE)
  flavornet = TRUE,          # 提取Flavornet数据 (默认: FALSE)
  synonyms = TRUE,           # 提取同义词 (默认: FALSE)
  uses = TRUE,               # 提取化合物用途 (默认: FALSE)
  verbose = TRUE             # 显示进度 (默认: TRUE)
)

# 使用ClassyFire添加化学分类
# 使用 ?extract_classyfire 查看详情
compounds_classified <- extract_classyfire(
  data = compounds_meta,     # 包含InChIKey列的数据框
  inchikey_col = "InChIKey", # InChIKey列名 (默认值)
  name_col = "Name"          # 用于进度显示的名称列 (默认值)
)
```

**关键参数说明:** - `name_col`, `cas_col`, `inchikey_col`:
指定包含标识符的列 (默认值: “Name”, “CAS”, “InChIKey”) - `verbose`: 设为
`FALSE` 可抑制进度信息 (默认: TRUE) - `cas`, `flavornet`, `synonyms`,
`uses`: 启用特定元数据类型的提取 (全部默认: FALSE) -
函数自动处理网络错误并提供详细的进度信息

</div>

------------------------------------------------------------------------

## 📄 **许可证**

本项目采用MIT许可证 - 详情请参阅 [LICENSE](LICENSE) 文件。

------------------------------------------------------------------------

## 👨‍🔬 **作者信息**

**苏启枝 (Qizhi Su)** - *包开发者* - 📧 邮箱: <sukissqz@gmail.com> - 🆔
ORCID: [0000-0002-8124-997X](https://orcid.org/0000-0002-8124-997X) - 🐙
GitHub: [@QizhiSu](https://github.com/QizhiSu)

------------------------------------------------------------------------

<div align="center">

**⭐ 如果您觉得labtools有用，请考虑给它一个星标！⭐**

[![GitHub
stars](https://img.shields.io/github/stars/QizhiSu/labtools?style=for-the-badge&logo=github&logoColor=white&labelColor=1a1a1a&color=ffd93d)](https://github.com/QizhiSu/labtools/stargazers)

------------------------------------------------------------------------

<h3 style="color: #00d4aa;">
🔬 让分析化学数据处理变得简单高效 🔬
</h3>
<p style="color: #a0a0a0; font-style: italic;">
Streamlining analytical chemistry workflows with precision and elegance
</p>

</div>
