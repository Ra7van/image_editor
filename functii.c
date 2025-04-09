#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
	double R, G, B;
} pixel;

void free_mem_matrice_pixel(pixel ***a, int linii)
{
	for (int i = 0; i < linii; i++) {
		free((*a)[i]);
		(*a)[i] = NULL;
	}
	free(*a);
	*a = NULL;
}

void free_mem_matrice_double(double ***a, int linii)
{
	for (int i = 0; i < linii; i++) {
		free((*a)[i]);
		(*a)[i] = NULL;
	}
	free(*a);
	*a = NULL;
}

void P2(FILE *in, double **a, int coloane, int linii)
{
	for (int i = 0; i < linii; i++)
		for (int j = 0; j < coloane; j++)
			fscanf(in, "%lf", &a[i][j]);
}

void P3(FILE *in, pixel **b, int coloane, int linii)
{
	for (int i = 0; i < linii; i++)
		for (int j = 0; j < coloane; j++) {
			fscanf(in, "%lf", &b[i][j].R);
			fscanf(in, "%lf", &b[i][j].G);
			fscanf(in, "%lf", &b[i][j].B);
		}
}

void P5(FILE *in, double **a, int coloane, int linii)
{
	unsigned char x;
	fread(&x, sizeof(unsigned char), 1, in);
	for (int i = 0; i < linii; i++) {
		for (int j = 0; j < coloane; j++) {
			fread(&x, sizeof(unsigned char), 1, in);
			a[i][j] = (double)x;
		}
	}
}

void P6(FILE *in, pixel **b, int coloane, int linii)
{
	unsigned char x;
	fread(&x, sizeof(unsigned char), 1, in);
	for (int i = 0; i < linii; i++) {
		for (int j = 0; j < coloane; j++) {
			fread(&x, sizeof(unsigned char), 1, in);
			b[i][j].R = (double)x;
			fread(&x, sizeof(unsigned char), 1, in);
			b[i][j].G = (double)x;
			fread(&x, sizeof(unsigned char), 1, in);
			b[i][j].B = (double)x;
		}
	}
}

int LOAD(FILE *in, FILE *inn, char type[], int coloane, int linii, double **a, pixel **b, int *ok, int *okp)
{
	switch (type[1]) {
	case '2':
		P2(inn, a, coloane, linii);
		(*ok) = 1;
		return 1;

	case '3':
		P3(inn, b, coloane, linii);
		(*okp) = 1;
		return 1;

	case '5':
		P5(in, a, coloane, linii);
		(*ok) = 1;
		return 1;

	case '6':
		P6(in, b, coloane, linii);
		(*okp) = 1;
		return 1;
	}
	return 0;
}

int SELECT(double **a, pixel **b, double **s1, pixel **s2, int *ok, int *okp, int x1, int y1, int x2, int y2)
{
	if ((*ok) == 1) {
		for (int i = x1; i < x2; i++)
			for (int j = y1; j < y2; j++)
				s1[j - y1][i - x1] = a[j][i];
	}
	if ((*okp) == 1) {
		for (int i = x1; i < x2; i++)
			for (int j = y1; j < y2; j++)
				s2[j - y1][i - x1] = b[j][i];
	}
	return 1;
}

void SELECT_ALL(double **a, pixel **b, double **s1, pixel **s2, int linii, int coloane, int *ok, int *okp)
{
	if ((*ok) == 1)
		for (int i = 0; i < linii; i++)
			for (int j = 0; j < coloane; j++)
				s1[i][j] = a[i][j];
	if ((*okp) == 1)
		for (int i = 0; i < linii; i++)
			for (int j = 0; j < coloane; j++)
				s2[i][j] = b[i][j];
}

void HISTOGRAM(double **a, double **s1, int linii, int coloane, int s1_linii, int s1_coloane, int x, int y)
{
	int nr;
	double stelute = 0, sum = 0, maxim = -1;
	int *v = (int *)calloc(256, sizeof(int));

	if (s1_linii == 0 || s1_coloane == 0) {
		for (int i = 0; i < linii; i++)
			for (int j = 0; j < coloane; j++)
				v[(int)(a[i][j])]++;
	} else {
		for (int i = 0; i < s1_linii; i++)
			for (int j = 0; j < s1_coloane; j++)
				v[(int)(s1[i][j])]++;
	}

	nr = 256 / y;
	int cnt = 256 / y;

	for (int i = 0; i < 256; i++) {
		if (i < nr)
			sum += v[i];
		if (i == nr - 1) {
			if (maxim < sum)
				maxim = sum;
			nr += cnt;
			sum = 0;
		}
	}

	sum = 0;
	nr = 256 / y;

	for (int i = 0; i < 256; i++) {
		if (i < nr)
			sum += v[i];
		if (i == nr - 1) {
			nr += cnt;
			stelute = sum / maxim * x;
			printf("%d\t|\t", (int)stelute);
			while ((int)stelute > 0) {
				printf("*");
				stelute--;
			}
			printf("\n");
			sum = 0;
		}
	}
	free(v);
	v = NULL;
	return;
}

void EQUALIZE(double **a, int linii, int coloane)
{
	int *f;
	double arie;
	f = (int *)calloc(257, sizeof(int));
	for (int i = 0; i < linii; i++)
		for (int j = 0; j < coloane; j++) {
			if (a[i][j] < 0)
				a[i][j] = 0;
			if (a[i][j] > 255)
				a[i][j] = 255;
			f[(int)a[i][j]]++;
		}
	for (int i = 1; i <= 256; i++)
		f[i] += f[i - 1];
	arie = linii * coloane;
	for (int i = 0; i < linii; i++)
		for (int j = 0; j < coloane; j++) {
			a[i][j] = 255 * (1.0 / arie) * f[(int)a[i][j]];
			if (a[i][j] < 0)
				a[i][j] = 0;
			if (a[i][j] > 255)
				a[i][j] = 255;
		}
	free(f);
	printf("Equalize done\n");
}

void rotateall_double(double **mat, double **a,  int *linii, int *coloane, int unghi)
{
	printf("Rotated %d\n", unghi);
	if (unghi == 90 || unghi == -270) {
		for (int i = 0; i < *linii; i++)
			for (int j = 0; j < *coloane; j++)
				a[j][*linii - 1 - i] = mat[i][j];
		free_mem_matrice_double(&mat, *linii);
		int aux = *linii;
		*linii = *coloane;
		*coloane = aux;
	} else if (unghi == 270 || unghi == -90) {
		for (int i = 0; i < *linii; i++)
			for (int j = 0; j < *coloane; j++)
				a[*coloane - j - 1][i] = mat[i][j];
		free_mem_matrice_double(&mat, *linii);
		int aux = *linii;
		*linii = *coloane;
		*coloane = aux;
	}
}

void rotate180_double(double **a, int linii, int coloane, int unghi)
{
	double **temp;
	temp = (double **)malloc(linii * sizeof(double *));
	for (int i = 0; i < linii; i++)
		temp[i] = (double *)malloc(coloane * sizeof(double));
	for (int i = 0; i < linii; i++)
		for (int j = 0; j < coloane; j++)
			temp[i][j] = a[i][j];
	for (int i = 0; i < linii; i++)
		for (int j = 0; j < coloane; j++)
			a[i][j] = temp[linii - 1 - i][coloane - 1 - j];
	free_mem_matrice_double(&temp, linii);
	printf("Rotated %d\n", unghi);
}

void rotateall_pixel(pixel **mat, pixel **b,  int *linii, int *coloane, int unghi)
{
	printf("Rotated %d\n", unghi);
	if (unghi == 90 || unghi == -270) {
		for (int i = 0; i < *linii; i++)
			for (int j = 0; j < *coloane; j++)
				b[j][*linii - 1 - i] = mat[i][j];
		free_mem_matrice_pixel(&mat, *linii);
		int aux = *linii;
		*linii = *coloane;
		*coloane = aux;
	} else if (unghi == 270 || unghi == -90) {
		for (int i = 0; i < *linii; i++)
			for (int j = 0; j < *coloane; j++)
				b[*coloane - j - 1][i] = mat[i][j];
		free_mem_matrice_pixel(&mat, *linii);
		int aux = *linii;
		*linii = *coloane;
		*coloane = aux;
	}
}

void rotate180_pixel(pixel **b, int linii, int coloane, int unghi)
{
	pixel **temp;
	temp = (pixel **)malloc(linii * sizeof(pixel *));
	for (int i = 0; i < linii; i++)
		temp[i] = (pixel *)malloc(coloane * sizeof(pixel));
	for (int i = 0; i < linii; i++)
		for (int j = 0; j < coloane; j++)
			temp[i][j] = b[i][j];
	for (int i = 0; i < linii; i++)
		for (int j = 0; j < coloane; j++)
			b[i][j] = temp[linii - 1 - i][coloane - 1 - j];
	free_mem_matrice_pixel(&temp, linii);
	printf("Rotated %d\n", unghi);
}

void rotate_double(double **a, int x1, int y1, int y2)
{
	double temp;
	int linie = y2 - y1;
	for (int i = 0; i < linie / 2; i++) {
		int start = i;
		int end = linie - i - 1;
		for (int j = start; j < end; j++) {
			int g = j - start;
			temp = a[start + y1][j + x1];
			a[start + y1][j + x1] = a[end + y1 - g][start + x1];
			a[end + y1 - g][start + x1] = a[end + y1][end + x1 - g];
			a[end + y1][end + x1 - g] = a[j + y1][end + x1];
			a[j + y1][end + x1] = temp;
		}
	}
}

void rotate_pixel(pixel **b, int x1, int y1, int y2)
{
	pixel temp;
	int linie = y2 - y1;
	for (int i = 0; i < linie / 2; i++) {
		int start = i;
		int end = linie - i - 1;
		for (int j = start; j < end; j++) {
			int g = j - start;
			temp = b[start + y1][j + x1];
			b[start + y1][j + x1] = b[end + y1 - g][start + x1];
			b[end + y1 - g][start + x1] = b[end + y1][end + x1 - g];
			b[end + y1][end + x1 - g] = b[j + y1][end + x1];
			b[j + y1][end + x1] = temp;
		}
	}
}

void rotate_d(double **a, int x1, int y1, int y2, int unghi)
{
	if (abs(unghi) == 360 || abs(unghi) == 0) {
		printf("Rotated %d\n", unghi);
		return;
	}
	if (unghi == 90 || unghi == -270)
		rotate_double(a, x1, y1, y2);
	if (unghi == 180 || unghi == -180)
		for (int i = 1; i <= 2; i++)
			rotate_double(a, x1, y1, y2);
	if (unghi == 270 || unghi == -90)
		for (int i = 1; i <= 3; i++)
			rotate_double(a, x1, y1, y2);
	printf("Rotated %d\n", unghi);
}

void rotate_p(pixel **b, int x1, int y1, int y2, int unghi)
{
	if (abs(unghi) == 360 || abs(unghi) == 0) {
		printf("Rotated %d\n", unghi);
		return;
	}
	if (unghi == 90 || unghi == -270)
		rotate_pixel(b, x1, y1, y2);
	if (unghi == 180 || unghi == -180)
		for (int i = 1; i <= 2; i++)
			rotate_pixel(b, x1, y1, y2);
	if (unghi == 270 || unghi == -90)
		for (int i = 1; i <= 3; i++)
			rotate_pixel(b, x1, y1, y2);
	printf("Rotated %d\n", unghi);
}

void CROP(double **a, pixel **b, double **s1, pixel **s2, int linii, int coloane, int *ok, int *okp)
{
	if ((*ok) == 1) {
		for (int i = 0; i < linii; i++)
			for (int j = 0; j < coloane; j++)
				a[i][j] = s1[i][j];
	}
	if ((*okp) == 1) {
		for (int i = 0; i < linii; i++)
			for (int j = 0; j < coloane; j++)
				b[i][j] = s2[i][j];
	}
}

int APPLY(pixel **b, int *okp, int x1, int x2, int y1, int y2, int linii, int coloane, char parametru[])
{
	double edge[3][3] = {{-1, -1, -1}, {-1, 8, -1}, {-1, -1, -1}};
	long sharpen[3][3] = {{0, -1, 0}, {-1, 5, -1}, {0, -1, 0}};
	double blur[3][3] = {{1.0 / 9, 1.0 / 9, 1.0 / 9},
	{1.0 / 9, 1.0 / 9, 1.0 / 9}, {1.0 / 9, 1.0 / 9, 1.0 / 9}};
	double gaussian_blur[3][3] = {{1.0 / 16, 2.0 / 16, 1.0 / 16},
	{2.0 / 16, 4.0 / 16, 2.0 / 16}, {1.0 / 16, 2.0 / 16, 1.0 / 16}};

	int valid = 0;
	double sum_R = 0, sum_G = 0, sum_B = 0;
	if ((*okp) != 1) {
		printf("Easy, Charlie Chaplin\n");
		return 0;
	}
	pixel **temp;
	temp = (pixel **)malloc(linii * sizeof(pixel *));
	for (int i = 0; i < linii; i++)
		temp[i] = (pixel *)malloc(coloane * sizeof(pixel));

	if (x1 == 0 && y1 == 0 && x2 == 0 && y2 == 0) {
		x1 = 1;
		y1 = 1;
		x2 = coloane - 1;
		y2 = linii - 1;
	} else {
		if (x1 == 0)
			x1 = 1;
		if (y1 == 0)
			y1 = 1;
		if (x2 == 0)
			x2 = 1;
		if (y2 == 0)
			y2 = 1;
		if (x2 == coloane)
			x2--;
		if (y2 == linii)
			y2--;
	}

	if (strcmp(parametru, "EDGE") == 0) {
		long sum_R = 0, sum_G = 0, sum_B = 0;
		for (int i = y1; i < y2; i++) {
			for (int j = x1; j < x2; j++) {
				sum_R = 0; sum_G = 0; sum_B = 0;
				for (int k = 0; k < 3; k++)
					for (int l = 0; l < 3; l++) {
						sum_R += b[i - 1 + k][j - 1 + l].R * edge[k][l];
						sum_G += b[i - 1 + k][j - 1 + l].G * edge[k][l];
						sum_B += b[i - 1 + k][j - 1 + l].B * edge[k][l];
					}
				if (sum_R < 0)
					sum_R = 0;
				else if (sum_R > 255)
					sum_R = 255;
				if (sum_G < 0)
					sum_G = 0;
				else if (sum_G > 255)
					sum_G = 255;
				if (sum_B < 0)
					sum_B = 0;
				else if (sum_B > 255)
					sum_B = 255;
				temp[i][j].R = sum_R;
				temp[i][j].G = sum_G;
				temp[i][j].B = sum_B;
			}
		}
		valid = 1;
	}

	if (strcmp(parametru, "SHARPEN") == 0) {
		sum_R = 0, sum_G = 0, sum_B = 0;
		for (int i = y1; i < y2; i++) {
			for (int j = x1; j < x2; j++) {
				sum_R = 0; sum_G = 0; sum_B = 0;
				for (int k = 0; k < 3; k++)
					for (int l = 0; l < 3; l++) {
						sum_R = sum_R + b[i - 1 + k][j - 1 + l].R * sharpen[k][l];
						sum_G = sum_G + b[i - 1 + k][j - 1 + l].G * sharpen[k][l];
						sum_B = sum_B + b[i - 1 + k][j - 1 + l].B * sharpen[k][l];
					}
				if (sum_R < 0)
					sum_R = 0;
				else if (sum_R > 255)
					sum_R = 255;
				if (sum_G < 0)
					sum_G = 0;
				else if (sum_G > 255)
					sum_G = 255;
				if (sum_B < 0)
					sum_B = 0;
				else if (sum_B > 255)
					sum_B = 255;
				temp[i][j].R = sum_R;
				temp[i][j].G = sum_G;
				temp[i][j].B = sum_B;
			}
		}
		valid = 1;
	}

	if (strcmp(parametru, "BLUR") == 0) {
		sum_R = 0, sum_G = 0, sum_B = 0;
		for (int i = y1; i < y2; i++) {
			for (int j = x1; j < x2; j++) {
				sum_R = 0; sum_G = 0; sum_B = 0;
				for (int k = 0; k < 3; k++)
					for (int l = 0; l < 3; l++) {
						sum_R += b[i - 1 + k][j - 1 + l].R * blur[k][l];
						sum_G += b[i - 1 + k][j - 1 + l].G * blur[k][l];
						sum_B += b[i - 1 + k][j - 1 + l].B * blur[k][l];
					}
				if (sum_R < 0)
					sum_R = 0;
				else if (sum_R > 255)
					sum_R = 255;
				if (sum_G < 0)
					sum_G = 0;
				else if (sum_G > 255)
					sum_G = 255;
				if (sum_B < 0)
					sum_B = 0;
				else if (sum_B > 255)
					sum_B = 255;
				temp[i][j].R = sum_R;
				temp[i][j].G = sum_G;
				temp[i][j].B = sum_B;
			}
		}
		valid = 1;
	}

	if (strcmp(parametru, "GAUSSIAN_BLUR") == 0) {
		sum_R = 0, sum_G = 0, sum_B = 0;
		for (int i = y1; i < y2; i++) {
			for (int j = x1; j < x2; j++) {
				sum_R = 0; sum_G = 0; sum_B = 0;
				for (int k = 0; k < 3; k++)
					for (int l = 0; l < 3; l++) {
						sum_R += b[i - 1 + k][j - 1 + l].R * gaussian_blur[k][l];
						sum_G += b[i - 1 + k][j - 1 + l].G * gaussian_blur[k][l];
						sum_B += b[i - 1 + k][j - 1 + l].B * gaussian_blur[k][l];
					}
				if (sum_R < 0)
					sum_R = 0;
				else if (sum_R > 255)
					sum_R = 255;
				if (sum_G < 0)
					sum_G = 0;
				else if (sum_G > 255)
					sum_G = 255;
				if (sum_B < 0)
					sum_B = 0;
				else if (sum_B > 255)
					sum_B = 255;
				temp[i][j].R = sum_R;
				temp[i][j].G = sum_G;
				temp[i][j].B = sum_B;
			}
		}
		valid = 1;
	}
	if (valid == 0) {
		printf("APPLY parameter invalid\n");
		free_mem_matrice_pixel(&temp, linii);
		return 0;
	} else {
		printf("APPLY %s done\n", parametru);
		for (int i = y1; i < y2; i++)
			for (int j = x1; j < x2; j++) {
				b[i][j].R = temp[i][j].R;
				b[i][j].G = temp[i][j].G;
				b[i][j].B = temp[i][j].B;
			}
		free_mem_matrice_pixel(&temp, linii);
		return 1;
	}
}

void SAVE(char name[], char type[], int linii, int coloane, int maxim, double **a, pixel **b, int *ok, int *okp)
{
	FILE *out = NULL;
	char *pozitie;
	if (strstr(name, "ascii")) {
		pozitie = strstr(name, " ascii");
		name[pozitie - name] = '\0';
		memmove(name, name + 1, strlen(name)); // pt a elimina spatiul de la inceput
		out = fopen(name, "wt");
		fprintf(out, "P");
		if (type[1] == '5')
			fprintf(out, "2\n");
		else if (type[1] == '6')
			fprintf(out, "3\n");
			else
				fprintf(out, "%d\n", type[1] - '0');
		fprintf(out, "%d %d\n", coloane, linii);
		fprintf(out, "%d\n", maxim);

		if ((*ok) == 1) {
			for (int i = 0; i < linii; i++) {
				for (int j = 0; j < coloane; j++)
					fprintf(out, "%ld ", (long)a[i][j]);
				fprintf(out, "\n");
			}
		}
		if ((*okp) == 1) {
			for (int i = 0; i < linii; i++) {
				for (int j = 0; j < coloane; j++) {
					fprintf(out, "%ld ", (long)b[i][j].R);
					fprintf(out, "%ld ", (long)b[i][j].G);
					fprintf(out, "%ld ", (long)b[i][j].B);
				}
				fprintf(out, "\n");
			}
		}
	} else {
		memmove(name, name + 1, strlen(name));
		out = fopen(name, "wb");
		fprintf(out, "P");
		if (type[1] == '2')
			fprintf(out, "5\n");
		else if (type[1] == '3')
			fprintf(out, "6\n");
			else
				fprintf(out, "%d\n", type[1] - '0');
		fprintf(out, "%d %d\n",  coloane, linii);
		fprintf(out, "%d\n", maxim);

		if ((*ok) == 1) {
			for (int i = 0; i < linii; i++) {
				for (int j = 0; j < coloane; j++) {
					unsigned char x = (int)a[i][j];
					fwrite(&x, sizeof(unsigned char), 1, out);
				}
			}
		}
		if ((*okp) == 1) {
			for (int i = 0; i < linii; i++) {
				for (int j = 0; j < coloane; j++) {
					unsigned char red = (int)b[i][j].R;
					unsigned char green = (int)b[i][j].G;
					unsigned char blue = (int)b[i][j].B;
					fwrite(&red, sizeof(unsigned char), 1, out);
					fwrite(&green, sizeof(unsigned char), 1, out);
					fwrite(&blue, sizeof(unsigned char), 1, out);
				}
			}
		}
	}
	fclose(out);
	printf("Saved %s\n", name);
}

void EXEC(void)
{
	FILE *in = NULL, *inn = NULL;
	double **a, **s1;
	pixel **b, **s2;
	char command[30], name[50], type[5], argument[30];
	int coloane = 0, linii = 0, maxim = -1, load = 0, *ok = malloc(sizeof(int)), *okp = malloc(sizeof(int));
	int s1_linii, s1_coloane, s2_linii, s2_coloane, citire = 0;
	int x1 = 0, x2 = 0, y1 = 0, y2 = 0, cx1, cx2, cy1, cy2;
	(*ok) = 0; (*okp) = 0;
	int crop = 0, select = 0;

	while (scanf("%s", command) == 1) {
		citire = 1;
		if (strcmp(command, "LOAD") == 0) {
			load = 0;
			if ((*ok) == 1)
				free_mem_matrice_double(&a, linii);
			if ((*okp) == 1)
				free_mem_matrice_pixel(&b, linii);
			if ((*ok) == 1 && s1_linii != 0)
				free_mem_matrice_double(&s1, s1_linii);
			if ((*okp) == 1 && s2_linii != 0)
				free_mem_matrice_pixel(&s2, s2_linii);
			(*ok) = 0; (*okp) = 0;
			scanf("%s", name);
			inn = fopen(name, "rt");
			if (!inn) {
				printf("Failed to load %s\n", name);
				load = 0;
			} else {
				fscanf(inn, "%s", type);
				fscanf(inn, "%d %d %d", &coloane, &linii, &maxim);
				x1 = 0; x2 = coloane; y1 = 0; y2 = linii;

				if (type[1] == '2' || type[1] == '5') {
					a = (double **)malloc(linii * sizeof(double *));
					for (int i = 0; i < linii; i++)
						a[i] = (double *)malloc(coloane * sizeof(double));
				}

				if (type[1] == '3' || type[1] == '6') {
					b = (pixel **)malloc(linii * sizeof(pixel *));
					for (int i = 0; i < linii; i++)
						b[i] = (pixel *)malloc(coloane * sizeof(pixel));
				}

				if (type[1] == '5' || type[1] == '6') {
					fclose(inn);
					inn = NULL;
					in = fopen(name, "rb");
					fscanf(in, "%s", type);
					fscanf(in, "%d %d %d", &coloane, &linii, &maxim);
				}

				load = LOAD(in, inn, type, coloane, linii, a, b, ok, okp);
				if (load == 0) {
					printf("Failed to load %s\n", name);
					continue;
				} else
					printf("Loaded %s\n", name);

				s1_linii = 0; s1_coloane = 0; s2_linii = 0; s2_coloane = 0;
				if (in)
					fclose(in);
				if (inn)
					fclose(inn);
				in = NULL; inn = NULL;
			}
			citire = 0;
		}
		if (strcmp(command, "SELECT") == 0) {
			int afisare = 0;
			fgets(argument, 30, stdin);
			memmove(argument, argument + 1, strlen(argument));
			argument[strlen(argument) - 1] = '\0';
			if (strstr(argument, "ALL") != 0) {
				if (load == 0) {
					printf("No image loaded\n");
					afisare = 1;
				} else {
					crop = 0; select = 1;
					x1 = 0; x2 = coloane; y1 = 0; y2 = linii;
					if ((*ok) == 1) {
						s1_linii = linii;
						s1_coloane = coloane;
						s1 = (double **)malloc(s1_linii * sizeof(double *));
						for (int i = 0; i < s1_linii; i++)
							s1[i] = (double *)malloc(s1_coloane * sizeof(double));
					}
					if ((*okp) == 1) {
						s2_linii = linii;
						s2_coloane = coloane;
						s2 = (pixel **)malloc(s2_linii * sizeof(pixel *));
						for (int i = 0; i < s2_linii; i++)
							s2[i] = (pixel *)malloc(s2_coloane * sizeof(pixel));
					}
					SELECT_ALL(a, b, s1, s2, linii, coloane, ok, okp);
					printf("Selected ALL\n");
					afisare = 1;
				}
				citire = 0;
			}
			if (sscanf(argument, "%d %d %d %d", &cx1, &cy1, &cx2, &cy2) == 4) {
				if (load == 0) {
					printf("No image loaded\n");
					afisare = 1;
				} else {
					int aux;
					if (cx1 > cx2) {
						aux = cx1;
						cx1 = cx2;
						cx2 = aux;
					}
					if (cy1 > cy2) {
						aux = cy1;
						cy1 = cy2;
						cy2 = aux;
					}
					if (cx1 < 0 || cx2 < 0 || cy1 < 0 || cy2 < 0 || cx1 > coloane || cx2 > coloane || cy1 > linii || cy2 > linii || (cx1 == cx2 || cy1 == cy2)) {
						printf("Invalid set of coordinates\n");
						afisare = 1;
					} else {
						select = 1; crop = 0;
						x1 = cx1; x2 = cx2; y1 = cy1; y2 = cy2;
						if ((*ok) == 1) {
							if (s1_linii != 0)
								free_mem_matrice_double(&s1, s1_linii);
							s1_linii = y2 - y1;
							s1_coloane = x2 - x1;

							s1 = (double **)malloc(s1_linii * sizeof(double *));
							for (int i = 0; i < s1_linii; i++)
								s1[i] = (double *)malloc(s1_coloane * sizeof(double));
						}

						if ((*okp) == 1) {
							if (s2_linii != 0)
								free_mem_matrice_pixel(&s2, s2_linii);
							s2_linii = y2 - y1;
							s2_coloane = x2 - x1;

							s2 = (pixel **)malloc(s2_linii * sizeof(pixel *));
							for (int i = 0; i < s2_linii; i++)
								s2[i] = (pixel *)malloc(s2_coloane * sizeof(pixel));
						}
						if (SELECT(a, b, s1, s2, ok, okp, x1, y1, x2, y2) == 1) {
							printf("Selected %d %d %d %d\n", x1, y1, x2, y2);
							afisare = 1;
						}
					}
				}
				citire = 0;
			}
			if (afisare == 0)
				printf("Invalid command\n");
		}
		if (strcmp(command, "HISTOGRAM") == 0) {
			char histo[20];
			int x, y, test;
			fgets(histo, sizeof(histo), stdin);
			memmove(histo, histo + 1, strlen(histo));
			histo[strlen(histo) - 1] = '\0';
			if (load == 0)
				printf("No image loaded\n");
			else {
				if (sscanf(histo, "%d %d %d", &x, &y, &test) != 2)
					printf("Invalid command\n");
				else if ((*okp) == 1 && (*ok) == 0)
					printf("Black and white image needed\n");
					else
						HISTOGRAM(a, s1, linii, coloane, s1_linii, s1_coloane, x, y);
			}
			citire = 0;
		}
		if (strcmp(command, "EQUALIZE") == 0) {
			if ((*okp) == 1 && (*ok) == 0)
				printf("Black and white image needed\n");
			else if ((*okp) == 0 && (*ok) == 0)
				printf("No image loaded\n");
			else
				EQUALIZE(a, linii, coloane);
			citire = 0;
		}
		if (strcmp(command, "ROTATE") == 0) {
			int unghi, test;
			char unghi_c[20];
			fgets(unghi_c, 20, stdin);
			memmove(unghi_c, unghi_c + 1, strlen(unghi_c));
			unghi_c[strlen(unghi_c) - 1] = '\0';
			if (load == 0)
				printf("No image loaded\n");
			else
				if (sscanf(unghi_c, "%d %d", &unghi, &test) != 1)
					printf("Invalid command\n");
				else
					if ((abs(unghi) % 90) != 0 || unghi < -360 || unghi > 360)
						printf("Unsupported rotation angle\n");
					else if (select == 1) {
						if ((x1 == 0 && y1 == 0 && x2 == coloane && y2 == linii) || (s1_linii == linii && s1_coloane == coloane)) {
							if ((*ok) == 1 && (*okp) == 0) {
								if (abs(unghi) == 90 || abs(unghi) == 270) {
									double **mat;
									mat = (double **)malloc(linii * sizeof(double *));
									for (int i = 0; i < linii; i++)
										mat[i] = (double *)malloc(coloane * sizeof(double));
									for (int i = 0; i < linii; i++)
										for (int j = 0; j < coloane; j++)
											mat[i][j] = a[i][j];
									free_mem_matrice_double(&a, linii);
									a = (double **)malloc(coloane * sizeof(double *));
									for (int i = 0; i < coloane; i++)
										a[i] = (double *)malloc(linii * sizeof(double));
									rotateall_double(mat, a, &linii, &coloane, unghi);
								} else if (abs(unghi) == 180)
									rotate180_double(a, linii, coloane, unghi);
									else
										if (abs(unghi) == 0 || abs(unghi) == 360)
											printf("Rotated %d\n", unghi);
							}
							if ((*ok) == 0 && (*okp) == 1) {
								if (abs(unghi) == 90 || abs(unghi) == 270) {
									pixel **mat;
									mat = (pixel **)malloc(linii * sizeof(pixel *));
									for (int i = 0; i < linii; i++)
										mat[i] = (pixel *)malloc(coloane * sizeof(pixel));
									for (int i = 0; i < linii; i++)
										for (int j = 0; j < coloane; j++)
											mat[i][j] = b[i][j];
									free_mem_matrice_pixel(&b, linii);
									b = (pixel **)malloc(coloane * sizeof(pixel *));
									for (int i = 0; i < coloane; i++)
										b[i] = (pixel *)malloc(linii * sizeof(pixel));
									rotateall_pixel(mat, b, &linii, &coloane, unghi);
								} else if (abs(unghi) == 180)
									rotate180_pixel(b, linii, coloane, unghi);
									else
										if (abs(unghi) == 0 || abs(unghi) == 360)
											printf("Rotated %d\n", unghi);
							}
						} else {
							if ((*ok) == 1 && (*okp) == 0) {
								if (s1_linii != s1_coloane && s1_linii != 0)
									printf("The selection must be square\n");
								else
									rotate_d(a, x1, y1, y2, unghi);
							}
							if ((*ok) == 0 && (*okp) == 1) {
								if (s2_linii != s2_coloane && s2_linii != 0)
									printf("The selection must be square\n");
								else
									rotate_p(b, x1, y1, y2, unghi);
							}
						}
					} else
						printf("Rotated %d\n", unghi);
			citire = 0;
		}
		if (strcmp(command, "CROP") == 0) {
			if (load == 0)
				printf("No image loaded\n");
			else if (crop == 0 && select == 1) {
				crop = 1; select = 0;
				if ((*ok) == 1)
					free_mem_matrice_double(&a, linii);
				if ((*okp) == 1)
					free_mem_matrice_pixel(&b, linii);
				if ((*ok) == 1 && s1_linii != 0) {
					linii = y2 - y1; coloane = x2 - x1;
					a = (double **)malloc(linii * sizeof(double *));
					for (int i = 0; i < linii; i++)
						a[i] = (double *)malloc(coloane * sizeof(double));

					s1_linii = 0; s1_coloane = 0;
				}
				if ((*okp) == 1 && s2_linii != 0) {
					linii = y2 - y1; coloane = x2 - x1;

					b = (pixel **)malloc(linii * sizeof(pixel *));
					for (int i = 0; i < linii; i++)
						b[i] = (pixel *)malloc(coloane * sizeof(pixel));
					s2_linii = 0; s2_coloane = 0;
				}
				CROP(a, b, s1, s2, linii, coloane, ok, okp);
				if ((*ok) == 1)
					free_mem_matrice_double(&s1, linii);
				if ((*okp) == 1)
					free_mem_matrice_pixel(&s2, linii);
				printf("Image cropped\n");
			} else
				printf("Image cropped\n");
			citire = 0;
		}
		if (strcmp(command, "APPLY") == 0) {
			char parametru[30];
			fgets(parametru, sizeof(parametru), stdin);
			if (load == 0)
				printf("No image loaded\n");
			else {
				memmove(parametru, parametru + 1, strlen(parametru));
				parametru[strlen(parametru) - 1] = '\0';
				if (parametru[0] == '\0')
					printf("Invalid command\n");
				else
					if (strcmp(parametru, "EDGE") != 0 && strcmp(parametru, "SHARPEN") != 0 && strcmp(parametru, "BLUR") != 0 && strcmp(parametru, "GAUSSIAN_BLUR") != 0)
						printf("APPLY parameter invalid\n");
					else
						if (APPLY(b, okp, x1, x2, y1, y2, linii, coloane, parametru) == 1) {
							free_mem_matrice_pixel(&s2, s2_linii);
							s2_linii = 0; s2_coloane = 0;
						}
			}
			citire = 0;
		}
		if (strcmp(command, "SAVE") == 0) {
			fgets(name, sizeof(name), stdin);
			name[strcspn(name, "\n")] = '\0';
			if (load == 0) {
				printf("No image loaded\n");
				if (in)
					fclose(in);
				if (inn)
					fclose(inn);
				in = NULL; inn = NULL;
			} else {
				SAVE(name, type, linii, coloane, maxim, a, b, ok, okp);
				if ((*okp) == 1 && s2_linii != 0) {
					free_mem_matrice_pixel(&s2, s2_linii);
					s2_linii = 0; s2_coloane = 0;
				}
				if ((*ok) == 1 && s1_linii != 0) {
					free_mem_matrice_double(&s1, s1_linii);
					s1_linii = 0; s1_coloane = 0;
				}
			}
			citire = 0;
		}
		if (strcmp(command, "EXIT") == 0) {
			if ((*ok) == 1 && linii != 0 && linii != s1_linii)
				free_mem_matrice_double(&a, linii);
			if ((*okp) == 1 && linii != 0 && linii != s2_linii)
				free_mem_matrice_pixel(&b, linii);
			if ((*ok) == 1 && s1_linii != 0)
				free_mem_matrice_double(&s1, s1_linii);
			if ((*okp) == 1 && s2_linii != 0)
				free_mem_matrice_pixel(&s2, s2_linii);
			if (in)
				fclose(in);
			if (inn)
				fclose(inn);
			in = NULL; inn = NULL;
			free(ok); free(okp);
			citire = 0;
			if (load == 0)
				printf("No image loaded\n");
			return;
		}
		if (citire == 1) {
			printf("Invalid command\n");
			int ch, o = 0;
			while ((ch = getchar()) != '\n' && ch != EOF)
				o++;
		}
	}
}
