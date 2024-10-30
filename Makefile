CUDA_PATH ?= /usr/local/cuda

NVCC = nvcc
NVCCFLAGS += -I$(CUDA_PATH)/include -L$(CUDA_PATH)/lib

SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = .

C_FILES = host.c \
	options.c \
	utils.c

CU_FILES = fdtd.cu device.cu

OBJ_FILES = $(C_FILES:.c=_c.o) $(CU_FILES:.cu=_cu.o)
OBJ_PATHS = $(addprefix $(OBJ_DIR)/, $(OBJ_FILES))

CC = gcc

GIT_COMMIT := $(shell git rev-parse HEAD | cut -c 1-8)

FLAGS += -I./inc -g
CCFLAGS +=  -DBUILD_VERSION="\"$(GIT_COMMIT)\""

ifneq (, $(shell which bear))
	BEAR := bear --append --
else
	BEAR :=
endif

$(shell mkdir -p $(OBJ_DIR) $(BIN_DIR))

all: fdtd

bear:
	$(BEAR) make -B fdtd

fdtd: $(OBJ_PATHS)
	echo $(OBJ_PATHS)
	$(NVCC) $(NVCCFLAGS) -o $@ $^

$(OBJ_DIR)/%_c.o: $(SRC_DIR)/%.c
	$(CC) $(CCFLAGS) $(FLAGS) -c $< -o $@

$(OBJ_DIR)/%_cu.o: $(SRC_DIR)/%.cu
	$(NVCC) $(NVCCFLAGS) $(FLAGS) -c $< -o $@

.PHONY: clean

clean:
	rm fdtd **/*.o

