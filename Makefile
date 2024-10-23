ifdef USE_HIP
NVCC = hipcc
NVCCFLAGS += --offload-arch=gfx90a
FLAGS += -D__HIP
else
NVCC = nvcc
NVCCFLAGS += -arch=sm_80 -I/usr/local/cuda/include -L/usr/local/cuda/lib
endif

SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = .

C_FILES = host.c \
	options.c \
	utils.c

CU_FILES = fdtd.cu

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

all:
	$(BEAR) make fdtd

.PHONY: clean fdtd time_detailed mem_check help

help:
	@echo "To compile for HIP, run `make USE_HIP=1`"

time_detailed: NVCCFLAGS+=-DTIME_DETAILED=1
time_detailed: fdtd

mem_check: NVCCFLAGS+=-DMEM_CHECK=1
mem_check: fdtd

fdtd: $(OBJ_PATHS)
	echo $(OBJ_PATHS)
	$(NVCC) $(NVCCFLAGS) -o $@ $^

$(OBJ_DIR)/%_c.o: $(SRC_DIR)/%.c
	$(CC) $(CCFLAGS) $(FLAGS) -c $< -o $@

$(OBJ_DIR)/%_cu.o: $(SRC_DIR)/%.cu FORCE
	$(NVCC) $(NVCCFLAGS) $(FLAGS) -c $< -o $@

FORCE:

clean:
	rm fdtd **/*.o


