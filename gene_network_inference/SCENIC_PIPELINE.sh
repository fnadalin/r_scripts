R_LIBS_USER="$HOME/.local/R-3.6.3":$R_LIBS_USER
export R_LIBS_USER

DB_DIR="$HOME/IIT/data/cisTopic_db"
export DB_DIR

WORKING_DIR=$( cd $( dirname $0 ) ; pwd )

# parameters
# OBJECT="$HOME/IIT/exp/analysis/perturb/object_filt_NEW.rds"
OBJECT="test/object.rds"
OUT_DIR="test"
SEED="123"
NUM_SEEDS="100"
NUM_CELLS="100"
CORES="1"

# extract the names of the cells
[[ -d ${OUT_DIR} ]] || mkdir ${OUT_DIR}
Rscript ${WORKING_DIR}/extract_cell_id_from_seurat_object.R \
	${OBJECT} \
	${OUT_DIR}/cell_id.txt

# generate random seeds
Rscript ${WORKING_DIR}/generate_random_seeds.R \
	${SEED} \
	${NUM_SEEDS} \
	${OUT_DIR}/${NUM_SEEDS}_random_seeds_seed${SEED}.txt

# generate the subset of cells
seed=$( head -1 ${OUT_DIR}/${NUM_SEEDS}_random_seeds_seed${SEED}.txt ) 
Rscript ${WORKING_DIR}/generate_random_cell_subset.R \
	${OUT_DIR}/cell_id.txt \
	${seed} \
	${NUM_CELLS} \
	${OUT_DIR}/cell_id_${NUM_CELLS}cells_seed${seed}.txt 

# generate the sparse sub-matrix
Rscript ${WORKING_DIR}/extract_sparse_submatrix_from_seurat_object.R \
	${OBJECT} \
	${OUT_DIR}/cell_id_${NUM_CELLS}cells_seed${seed}.txt \
	${OUT_DIR}/${NUM_CELLS}cells_seed${seed}

# SCENIC - create regulons on the sub-matrix
Rscript ${WORKING_DIR}/scenic_create_regulons.R \
	${OUT_DIR}/${NUM_CELLS}cells_seed${seed} \
	${OUT_DIR}/${NUM_CELLS}cells_seed${seed}/scenic \
	${CORES} 



exit



