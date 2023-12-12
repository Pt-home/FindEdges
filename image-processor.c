#include <emscripten.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef unsigned char byte;
typedef unsigned int uint32;
typedef unsigned long long uint64;
typedef float float32;

// Функция для клампинга значений
float clamp(float x, float lower, float upper) {
    if (x < lower) return lower;
    if (x > upper) return upper;
    return x;
}

#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef unsigned char byte;

typedef struct {
    int** values;
    int size;
} Kernel;

Kernel* createKernel(int size) {
    Kernel* k = (Kernel*) malloc(sizeof(Kernel));
    k->values = (int**) malloc(size * sizeof(int*));
    for (int i = 0; i < size; ++i) {
        k->values[i] = (int*) malloc(size * sizeof(int));
    }
    k->size = size;
    return k;
}

void freeKernel(Kernel* k) {
    for (int i = 0; i < k->size; ++i) {
        free(k->values[i]);
    }
    free(k->values);
    free(k);
}

void applyFE(byte* src, int width, int height, const char* method) {
    Kernel* Gx;
    Kernel* Gy;

    if (strcmp(method, "scharr") == 0) {        
        Gx = createKernel(3);
        Gy = createKernel(3);

        int tGx[3][3] = { { -3, 0, 3 }, { -10, 0, 10 }, { -3, 0, 3 } };
        int tGy[3][3] = { { -3, -10, -3 }, { 0, 0, 0 }, { 3, 10, 3 } };
        for (int i = 0; i < 3; ++i) {
            memcpy(Gx->values[i], tGx[i], 3 * sizeof(int));
            memcpy(Gy->values[i], tGy[i], 3 * sizeof(int));
        }
    } 
    else if (strcmp(method, "sobel") == 0) {        
        Gx = createKernel(3);
        Gy = createKernel(3);

        int tGx[3][3] = { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } };
        int tGy[3][3] = { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } };
        for (int i = 0; i < 3; ++i) {
            memcpy(Gx->values[i], tGx[i], 3 * sizeof(int));
            memcpy(Gy->values[i], tGy[i], 3 * sizeof(int));
        }
    }
    else if (strcmp(method, "prewitt") == 0) {        
        Gx = createKernel(3);
        Gy = createKernel(3);

        int tGx[3][3] = { { -1, 0, 1 }, { -1, 0, 1 }, { -1, 0, 1 } };
        int tGy[3][3] = { { -1, -1, -1 }, { 0, 0, 0 }, { 1, 1, 1 } };
        for (int i = 0; i < 3; ++i) {
            memcpy(Gx->values[i], tGx[i], 3 * sizeof(int));
            memcpy(Gy->values[i], tGy[i], 3 * sizeof(int));
        }
    }
    else if (strcmp(method, "roberts") == 0) { 
        Gx = createKernel(2);
        Gy = createKernel(2);

        int tGx[2][2] = { { 0, 1 }, { -1, 0 } };
        int tGy[2][2] = { { 1, 0 }, { 0, -1 } };
        for (int i = 0; i < 2; ++i) {
            memcpy(Gx->values[i], tGx[i], 2 * sizeof(int));
            memcpy(Gy->values[i], tGy[i], 2 * sizeof(int));
        }
    } 

    byte* edges = (byte*)malloc(width * height * 4);
    int* tempGx = (int*)calloc(width * height * 3, sizeof(int));
    int* tempGy = (int*)calloc(width * height * 3, sizeof(int));

    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            int sumX[3] = {0, 0, 0};
            int sumY[3] = {0, 0, 0};

            for (int i = -1; i <= Gx->size - 2; i++) {
                for (int j = -1; j <= Gx->size - 2; j++) {
                    byte* curPixel = src + ((y + i) * width + (x + j)) * 4;
                    for (int k = 0; k < 3; k++) {
                        sumX[k] += curPixel[k] * Gx->values[i + 1][j + 1];
                        sumY[k] += curPixel[k] * Gy->values[i + 1][j + 1];
                    }
                }
            }

            for (int k = 0; k < 3; k++) {
                tempGx[(y * width + x) * 3 + k] = sumX[k];
                tempGy[(y * width + x) * 3 + k] = sumY[k];
            }
        }
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            byte* curPixel = edges + (y * width + x) * 4;
            for (int k = 0; k < 3; k++) {
                int val = sqrt(tempGx[(y * width + x) * 3 + k] * tempGx[(y * width + x) * 3 + k] + 
                               tempGy[(y * width + x) * 3 + k] * tempGy[(y * width + x) * 3 + k]);
                if (val > 255) val = 255;
                curPixel[k] = (byte)val;
            }
            curPixel[3] = 255;
        }
    }

    memcpy(src, edges, width * height * 4);

    free(edges);
    free(tempGx);
    free(tempGy);
    freeKernel(Gx);
    freeKernel(Gy);
}

void applyLaplacian(byte* src, int width, int height, int kernelSize) {
    
    int size = kernelSize * kernelSize;
    int* laplacianKernel = (int*)malloc(size * sizeof(int));
    
    if (kernelSize == 3) {
        // Оператор Лапласа 3x3
        int tKernel[9] = { 0, 1, 0, 1, -4, 1, 0, 1, 0 };
        memcpy(laplacianKernel, tKernel, sizeof(tKernel));
    } else if (kernelSize == 5) {
        // Оператор Лапласа 5x5
        int tKernel[25] = { 0, 0, 1, 0, 0, 0, 1, 2, 1, 0, 1, 2, -16, 2, 1, 0, 1, 2, 1, 0, 0, 0, 1, 0, 0 };
        memcpy(laplacianKernel, tKernel, sizeof(tKernel));
    }

    byte* edges = (byte*)malloc(width * height * 4);
    int* temp = (int*)calloc(width * height * 3, sizeof(int));
    
    int offset = kernelSize / 2;

    for (int y = offset; y < height - offset; y++) {
        for (int x = offset; x < width - offset; x++) {
            int sum[3] = {0, 0, 0};

            // Apply Laplacian operator to each pixel separately, handling each RGB channel separately
            for (int i = -offset; i <= offset; i++) {
                for (int j = -offset; j <= offset; j++) {
                    byte* curPixel = src + ((y + i) * width + (x + j)) * 4;
                    int curKernelElem = laplacianKernel[(i + offset) * kernelSize + (j + offset)];

                    for (int c = 0; c < 3; c++) {
                        sum[c] += curPixel[c] * curKernelElem;
                    }
                }
            }

            // Save the results
            int tempIndex = (y * width + x) * 3;
            for (int c = 0; c < 3; c++) {
                temp[tempIndex + c] = sum[c];
            }
        }
    }

    // Copy the results to the edges buffer
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int tempIndex = (y * width + x) * 3;
            int edgeIndex = (y * width + x) * 4;

            for (int c = 0; c < 3; c++) {
                int value = temp[tempIndex + c];

                // Clamping the value to the range [0, 255]
                if (value < 0)
                    value = 0;
                else if (value > 255)
                    value = 255;

                edges[edgeIndex + c] = (byte)value;
            }

            // Copy alpha channel from the source image
            edges[edgeIndex + 3] = src[edgeIndex + 3];
        }
    }
    
    memcpy(src, edges, height * width * 4);

    // Free the allocated buffers
    free(laplacianKernel);
    free(temp);
    free(edges);
}

void applyKirsch(byte* src, int width, int height) {
    int kirschKernel[8][3][3] = {
        { {-3, -3, 5}, {-3, 0, 5}, {-3, -3, 5} },
        { {-3, 5, 5}, {-3, 0, 5}, {-3, -3, -3} },
        { {5, 5, 5}, {-3, 0, -3}, {-3, -3, -3} },
        { {5, 5, -3}, {5, 0, -3}, {-3, -3, -3} },
        { {5, -3, -3}, {5, 0, -3}, {5, -3, -3} },
        { {-3, -3, -3}, {5, 0, -3}, {5, 5, -3} },
        { {-3, -3, -3}, {-3, 0, -3}, {5, 5, 5} },
        { {-3, -3, -3}, {-3, 0, 5}, {-3, 5, 5} }
    };

    byte* edges = (byte*)malloc(width * height * 4);
    int* temp = (int*)calloc(width * height * 3, sizeof(int));
    
    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            int sum[8][3] = {0};

            // Apply Kirsch operator to each pixel separately, handling each RGB channel separately
            for (int k = 0; k < 8; k++) {
                for (int i = -1; i <= 1; i++) {
                    for (int j = -1; j <= 1; j++) {
                        byte* curPixel = src + ((y + i) * width + (x + j)) * 4;
                        int curKernelElem = kirschKernel[k][i+1][j+1];

                        for (int c = 0; c < 3; c++) {
                            sum[k][c] += curPixel[c] * curKernelElem;
                        }
                    }
                }
            }

            // Save the maximum of the results
            int tempIndex = (y * width + x) * 3;
            for (int c = 0; c < 3; c++) {
                int max = sum[0][c];
                for (int k = 1; k < 8; k++) {
                    if (sum[k][c] > max)
                        max = sum[k][c];
                }
                temp[tempIndex + c] = max;
            }
        }
    }

    // Copy the results to the edges buffer
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int tempIndex = (y * width + x) * 3;
            int edgeIndex = (y * width + x) * 4;

            for (int c = 0; c < 3; c++) {
                int value = temp[tempIndex + c];

                // Clamping the value to the range [0, 255]
                if (value < 0)
                    value = 0;
                else if (value > 255)
                    value = 255;

                edges[edgeIndex + c] = (byte)value;
            }

            // Copy alpha channel from the source image
            edges[edgeIndex + 3] = src[edgeIndex + 3];
        }
    }
    
    memcpy(src, edges, height * width * 4);

    // Free the allocated buffers
    free(temp);
    free(edges);
}

void applyRobinson(byte* src, int width, int height) {
    int G[8][3][3] = {
        { { -1, 0, 1 }, { -2, 0, 2 }, { -1, 0, 1 } },   // N
        { { 0, 1, 2 }, { -1, 0, 1 }, { -2, -1, 0 } },   // NE
        { { 1, 2, 1 }, { 0, 0, 0 }, { -1, -2, -1 } },   // E
        { { 2, 1, 0 }, { 1, 0, -1 }, { 0, -1, -2 } },   // SE
        { { 1, 0, -1 }, { 2, 0, -2 }, { 1, 0, -1 } },   // S
        { { 0, -1, -2 }, { 1, 0, -1 }, { 2, 1, 0 } },   // SW
        { { -1, -2, -1 }, { 0, 0, 0 }, { 1, 2, 1 } },   // W
        { { -2, -1, 0 }, { -1, 0, 1 }, { 0, 1, 2 } }    // NW
    };

    byte* edges = (byte*)malloc(width * height * 4);
    int* tempG[8];
    for (int i = 0; i < 8; i++) {
        tempG[i] = (int*)calloc(width * height * 3, sizeof(int));
    }

    for (int y = 1; y < height - 1; y++) {
        for (int x = 1; x < width - 1; x++) {
            int sumG[8][3] = {0};

            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    byte* curPixel = src + ((y + i) * width + (x + j)) * 4;

                    for (int g = 0; g < 8; g++) {
                        for (int c = 0; c < 3; c++) {
                            sumG[g][c] += curPixel[c] * G[g][i+1][j+1];
                        }
                    }
                }
            }

            int tempIndex = (y * width + x) * 3;
            for (int g = 0; g < 8; g++) {
                for (int c = 0; c < 3; c++) {
                    tempG[g][tempIndex + c] = sumG[g][c];
                }
            }
        }
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int edgeIndex = (y * width + x) * 4;
            int maxMagnitude[3] = {0};

            for (int g = 0; g < 8; g++) {
                int tempIndex = (y * width + x) * 3;

                for (int c = 0; c < 3; c++) {
                    int magnitude = abs(tempG[g][tempIndex + c]);
                    if (magnitude > maxMagnitude[c]) {
                        maxMagnitude[c] = magnitude;
                    }
                }
            }

            for (int c = 0; c < 3; c++) {
                edges[edgeIndex + c] = (byte)maxMagnitude[c];
            }

            edges[edgeIndex + 3] = src[edgeIndex + 3];
        }
    }

    memcpy(src, edges, height * width * 4);

    for (int i = 0; i < 8; i++) {
        free(tempG[i]);
    }
    free(edges);
}

EMSCRIPTEN_KEEPALIVE
byte* wasmAlloc(uint32 width, uint32 height) {
  return malloc(width * height * 4);
}

EMSCRIPTEN_KEEPALIVE
void wasmFree(byte* p) {
  free(p);
}

EMSCRIPTEN_KEEPALIVE
void wasmProcess(byte* inputData, uint32 size, int width, int height, const char* method ) {
	
	if (strcmp(method, "laplace3") == 0) {
		applyLaplacian(inputData, width, height, 3);
	} else if (strcmp(method, "laplace5") == 0) {
		applyLaplacian(inputData, width, height, 5);
	} else if (strcmp(method, "kirsch") == 0) {
		applyKirsch(inputData, width, height);
	} else if (strcmp(method, "robinson") == 0) {
		applyRobinson(inputData, width, height);
	} else {
		applyFE(inputData, width, height, method);
	}
	
}



