typedef struct {
    double R, G, B;
} pixel;

void free_mem_matrice_pixel(pixel ***a, int linii);
void free_mem_matrice_double(double ***a, int linii);
void P2(FILE *in, double **a, int coloane, int linii);
void P3(FILE *in, pixel **b, int coloane, int linii);
void P5(FILE *in, double **a, int coloane, int linii);
void P6(FILE *in, pixel **b, int coloane, int linii);
int LOAD (FILE *in, FILE *inn, char type[], int coloane, int linii, double **a, pixel **b, int *ok, int *okp);
int SELECT(double **a, pixel **b, double **s1, pixel **s2, int *ok, int *okp, int x1, int y1, int x2, int y2);
void SELECT_ALL(double **a, pixel **b, double **s1, pixel **s2, int linii, int coloane, int *ok, int *okp);
void HISTOGRAM(double **a, double **s1, int linii, int coloane, int s1_linii, int s1_coloane, int x, int y);
void EQUALIZE(double **a, int linii, int coloane);
void rotateall_double(double **mat, double **a,  int *linii, int *coloane, int unghi);
void rotate180_double(double **a, int linii, int coloane, int unghi);
void rotateall_pixel(pixel **mat, pixel **b,  int *linii, int *coloane, int unghi);
void rotate180_pixel(pixel **b, int linii, int coloane, int unghi);
void rotate_double(double **a, int x1, int y1, int y2);
void rotate_pixel(pixel **b, int x1, int y1, int y2);
void rotate_d(double **a, int x1, int y1, int y2, int unghi);
void rotate_p(pixel **b, int x1, int y1, int y2, int unghi);
void CROP(double **a, pixel **b, double **s1, pixel **s2, int linii, int coloane, int *ok, int *okp);
int APPLY(pixel **b, int *okp, int x1, int x2, int y1, int y2, int linii, int coloane, char parametru[]);
void SAVE(char name[], char type[], int linii, int coloane, int maxim, double **a, pixel **b, int *ok, int *okp);
void EXEC();
