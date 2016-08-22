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
#define OVERLAP 1

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

// double maxAmp, minAmp;
int q,i;

float freq2rate (int frequency)
{
  return (float)frequency/((float)SAMPLING_RATE/2);
}

void reduceFreq(fftw_complex *spectrum, int from, int to, float rate)
{
  for (i = from; i < to; i += 1)
  {
    spectrum[i] *= rate;
  }
}

void reducePitch(fftw_complex* spectrum, int shift)
{
  for (i = 0; i < SPECTRUM_SIZE; i++)
  {
    if(i >= SPECTRUM_SIZE - shift)
    {
      spectrum[i] = 0;
    }
    else
    {
      spectrum[i] = spectrum[i+shift];
    }
  }
}

long test1;
long preOver, over;

int main(int argc, char **argv){

  //Inputs
  float pre_amp = strtof(argv[3],NULL);
  float reduce_low_rate = strtof(argv[4],NULL);
  float cut_low_rate = freq2rate(atoi(argv[5]));
  float reduce_mid_rate = strtof(argv[6],NULL);
  float cut_high_rate = freq2rate(atoi(argv[7]));
  float reduce_high_rate = strtof(argv[8],NULL);
  int pitch_shift = atoi(argv[9]);
  int shift_val = atoi(argv[10]);
  float shift_mul = strtof(argv[11],NULL);
  float post_amp = strtof(argv[12],NULL);

  //Inputs check
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
  // long long *out_buf=malloc(10000000);
  int16_t *out_buf=malloc(10000000);
  int size=fread(buf,1,10000000,f)/2;

  fftw_complex *mid = fftw_alloc_complex(WINDOW);
  double *in = malloc(WINDOW*sizeof(double));
  double *out = malloc(WINDOW*sizeof(double));
  fftw_plan planForward = fftw_plan_dft_r2c_1d(WINDOW, in, mid, FFTW_ESTIMATE);
  fftw_plan planBackward = fftw_plan_dft_c2r_1d(WINDOW, mid, out, FFTW_ESTIMATE);
  // fftw_plan planBackward = fftw_plan_dft_c2r_1d(WINDOW, mid1, out, FFTW_ESTIMATE);

  //Pre-amp
  int maxAmp = 0;
  for (i = 44; i < size; i++)
  {
    buf[i] *= pre_amp;
    if(abs(buf[i]) > maxAmp)
    {
      maxAmp = abs(buf[i]);
    }
  }

  int t = 0;
  for (q = 44; q + WINDOW * (OVERLAP) <= size; q += WINDOW * (OVERLAP))
  {
    for (i = 0; i < WINDOW; i++)
    {
      // in[i] = (double)buf[q+i] / MAX_INT * 0.5 * (1 - cos(2*PI*i/(WINDOW-1)));//Hann function
      // in[i] = (double)buf[q+i] * 0.5 * (1 - cos(2*PI*i/(WINDOW-1)));//Hann function
      in[i] = (double)buf[q+i];//Rect function
    }

    fftw_execute(planForward);

    reduceFreq(mid,0,(int)(cut_low_rate*SPECTRUM_SIZE),reduce_low_rate);
    reduceFreq(mid,(int)(cut_low_rate*SPECTRUM_SIZE),(int)(cut_high_rate*SPECTRUM_SIZE),reduce_mid_rate);
    reduceFreq(mid,(int)(cut_high_rate*SPECTRUM_SIZE),(SPECTRUM_SIZE),reduce_high_rate);

    reducePitch(mid, pitch_shift);

    fftw_execute(planBackward);

    long test;
    for (i = 0; i < WINDOW; i++)
    {
      test = out_buf[q+i];
      out_buf[q+i] += (out[i] / WINDOW * OVERLAP);
      if(test < 0 && out_buf[q+i] > 0)
      {
        printf("%s\n", "FUCK YOU");
        test1++;
      }
      else if (test > 0 && out_buf[q+i] < 0)
      {
        printf("%s\n", "VERY MUCH");
        test1++;
      }
      if((out[i] / WINDOW * OVERLAP) > MAX_INT * OVERLAP)
      {
        preOver++;
      }
      if((out[i] / WINDOW * OVERLAP) > MAX_INT)
      {
        over++;
      }
    }
  }
  printf("%s %ld\n", "ERROR", test1);
  printf("%s %ld\n", "ERROR preOver", preOver);
  printf("%s %ld\n", "ERROR over", over);

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
    out_buf[i] *= 0.5;
    out_buf[i] += buf[i - shift_val] * shift_mul * 0.5;
  }

  int postMaxAmp = 0;
  for (i = 44; i < size; i++)
  {
    if (labs(out_buf[i]) > postMaxAmp)
    {
      postMaxAmp = labs(out_buf[i]);
    }
  }

  double normMul = (double)maxAmp/postMaxAmp;
  for (i = 44; i < size; i++)
  {
    out_buf[i] *= normMul;
    buf[i] = (int16_t)out_buf[i];
  }

  fftw_destroy_plan(planForward);
  fftw_destroy_plan(planBackward);
  free(in);
  fftw_free(mid);
  free(out);

  FILE *fo=fopen(argv[2],"wb");
  fwrite(buf,2,size,fo);
  fclose(fo);

  return 0;
}
