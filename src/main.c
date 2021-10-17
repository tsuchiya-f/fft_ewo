#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "HF_EWO_fft.h"

#define N (256)
#define NW (34)

float f_wave[N * 2];
float f_spec[N * 2];
int i_workarea[NW];
float f_sin_table[N];

float xin[N];
float yin[N];

void main ( void ){

	int i;
	float f_rad;
	float pow, pha, rel, img;

	FILE *fp_w, *fp_s;

	fp_w = fopen("waveform.dat", "w");
	fp_s = fopen("spectrum.dat", "w");

	for (i = 0; i < N; i++)
	{
		f_rad = (float)i / 8.0 * 2.0 * 3.14159265;
		//xin[i] = sin( f_rad );
		//yin[i] = 0.0;

		//xin[i] = cos(f_rad);
		//yin[i] = 0.0;

		//xin[i] = 0.0;
		//yin[i] = sin(f_rad);

		xin[i] = 0.0;
		yin[i] = cos(f_rad);

		//		xin[i] = sin(f_rad);
		//		yin[i] = sin(f_rad + 3.14159265 * 0.5);
	}

	// HF_EWO_Hanning ( xin, xin, N);
	// HF_EWO_Hanning ( yin, yin, N);
	// HF_EWO_Hamming(xin, xin, N);
	// HF_EWO_Hamming(yin, yin, N);
	// HF_EWO_Blackman(xin, N);
	// HF_EWO_Blackman(yin, N);

	for (i = 0; i < N; i++)
	{
		f_wave[i * 2] = xin[i];
		f_wave[i * 2 + 1] = yin[i];
	}

	i_workarea[0] = 0;
	HF_EWO_cdft(N * 2, 1, f_wave, i_workarea, f_sin_table);

	memcpy(f_spec, f_wave, N*2*sizeof(float));

	i_workarea[0] = 0;
	HF_EWO_cdft(N * 2, -1, f_wave, i_workarea, f_sin_table);

	fprintf(fp_s, "# real img, power, phase\n");
	for (i = 0; i < N; i++)
	{
		rel = f_spec[i * 2];
		img = f_spec[i * 2 + 1];
		pow = sqrt(rel*rel + img*img);
		pha = atan2(img, rel) * 180.0 / 3.14159265;
		fprintf ( fp_s, "%d %f %f %f %f\n", i, rel, img, pow, pha );
	}
	fclose(fp_s);

	fprintf(fp_w, "# x(org)  y(org)  x(ifft)  y(ifft)\n");
	for (i = 0; i < N; i++)
	{
		fprintf(fp_w, "%d %f %f %f %f \n", i, xin[i], yin[i], f_wave[i * 2], f_wave[i * 2 + 1]);
	}
	fclose(fp_w);

}
