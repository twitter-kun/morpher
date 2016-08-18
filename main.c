#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define REAL 0
#define IMAG 1

#define WINDOW 2048
#define SPECTRUM_SIZE WINDOW/2 + 1
#define SAMPLING_RATE 44100
#define OVERLAP 0.1

// #define CUT_LOW_RATE 0.03
// #define REDUCE_RATE 0.04
//
// #define SHIFT_VAL 500
// #define SHIFT_MUL 0.6
//
// #define PRE_AMP 0.2
// #define AMP 8

#define PI acos(-1.0)
#define MAX_INT 32767

double maxAmp, minAmp;

float freq2rate (int frequency)
{
  return (float)frequency/((float)SAMPLING_RATE/2);
}

int main(int argc, char **argv){

  float pre_amp = strtof(argv[3],NULL);
  float reduce_low_rate = strtof(argv[4],NULL);
  float cut_low_rate = freq2rate(atoi(argv[5]));
  float reduce_mid_rate = strtof(argv[6],NULL);
  float cut_high_rate = freq2rate(atoi(argv[7]));
  float reduce_high_rate = strtof(argv[8],NULL);
  int shift_val = atoi(argv[9]);
  float shift_mul = strtof(argv[10],NULL);
  float post_amp = strtof(argv[11],NULL);

  if(cut_high_rate < cut_low_rate ||
     reduce_high_rate < 0 ||
     reduce_mid_rate < 0 ||
     reduce_low_rate < 0)
  {
    printf("INPUT ERROR\n");
    return 1;
  }

  FILE *f=fopen(argv[1],"rb");
  int16_t *buf=malloc(10000000);
  int16_t *out_buf=malloc(10000000);
  int size=fread(buf,1,10000000,f)/2;

  fftw_complex *mid = fftw_alloc_complex(WINDOW);
  double *in = malloc(WINDOW*sizeof(double));
  double *out = malloc(WINDOW*sizeof(double));
  fftw_plan planForward = fftw_plan_dft_r2c_1d(WINDOW, in, mid, FFTW_ESTIMATE);
  fftw_plan planBackward = fftw_plan_dft_c2r_1d(WINDOW, mid, out, FFTW_ESTIMATE);

  int q,i;
  for (i = 44; i < size; i++)
  {
    buf[i] *= pre_amp;
  }

  int t = 0;
  for (q = 44; q + WINDOW * (OVERLAP) <= size; q += WINDOW * (OVERLAP))
  {
    for (i = 0; i < WINDOW; i++)
    {
      // in[i] = (double)buf[q+i] / MAX_INT * 0.5 * (1 - cos(2*PI*i/(WINDOW-1)));//Hann function
      in[i] = (double)buf[q+i] * 0.5 * (1 - cos(2*PI*i/(WINDOW-1)));//Hann function
      // in[i] = (double)buf[q+i];//Rect function
    }

    fftw_execute(planForward);
    // if (q >= 21400 && t != 1)
    // {
    //   t = 1;
    //   for (i = 0; i < WINDOW; i++)
    //   {
    //     printf("%f\t", cabs(mid[i]));
    //   }
    //   printf("\n");
    // }

    // for (i = SPECTRUM_SIZE; i < WINDOW; i++)
    // {
    //   if(cabs(mid[i]) < minAmp)
    //     minAmp = cabs(mid[i]);
    //   else if (cabs(mid[i]) > maxAmp)
    //     maxAmp = cabs(mid[i]);
    // }

    for (i = 0; i < (int)(cut_low_rate*SPECTRUM_SIZE); i += 1)
    {
      mid[i] *= reduce_low_rate;
    }

    for (i = (int)(cut_low_rate*SPECTRUM_SIZE); i < (int)(cut_high_rate*SPECTRUM_SIZE); i += 1)
    {
      mid[i] *= reduce_mid_rate;
    }

    for (i = (int)(cut_high_rate*SPECTRUM_SIZE); i < (SPECTRUM_SIZE); i += 1)
    {
      mid[i] *= reduce_high_rate;
    }

    fftw_execute(planBackward);

    for (i = 0; i < WINDOW; i++)
    {
      out_buf[q+i] += (int16_t)(out[i] / WINDOW * OVERLAP);
    }
  }

  for (i = 0; i < 44; i++)
  {
    out_buf[i] = buf[i];
  }

  for (i = 44; i < size; i++)
  {
    out_buf[i] *= post_amp;
  }

  for (i = 44 + shift_val; i < size; i++)
  {
    out_buf[i] += buf[i - shift_val] * shift_mul;
    out_buf[i] *= 0.5;
  }

  fftw_destroy_plan(planForward);
  fftw_destroy_plan(planBackward);
  free(in);
  fftw_free(mid);
  free(out);

  FILE *fo=fopen(argv[2],"wb");
  fwrite(out_buf,2,size,fo);
  fclose(fo);

  printf("Min: %f, Max: %f\n",minAmp,maxAmp);

  return 0;
}
