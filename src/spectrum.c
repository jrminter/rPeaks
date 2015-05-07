// Author: Miroslav Morhac   27/05/99

//__________________________________________________________________________
//   THIS CLASS CONTAINS ADVANCED SPECTRA PROCESSING FUNCTIONS.            //
//                                                                         //
//   ONE-DIMENSIONAL BACKGROUND ESTIMATION FUNCTIONS                       //
//   ONE-DIMENSIONAL SMOOTHING FUNCTIONS                                   //
//   ONE-DIMENSIONAL DECONVOLUTION FUNCTIONS                               //
//   ONE-DIMENSIONAL PEAK SEARCH FUNCTIONS                                 //
//                                                                         //
//   These functions were written by:                                      //
//   Miroslav Morhac                                                       //
//   Institute of Physics                                                  //
//   Slovak Academy of Sciences                                            //
//   Dubravska cesta 9, 842 28 BRATISLAVA                                  //
//   SLOVAKIA                                                              //
//                                                                         //
//   email:fyzimiro@savba.sk,    fax:+421 7 54772479                       //
//                                                                         //
//  The original code in C has been repackaged as a C++ class by R.Brun    //
//                                                                         //
//  The algorithms in this class have been published in the following      //
//  references:                                                            //
//   [1]  M.Morhac et al.: Background elimination methods for              //
//   multidimensional coincidence gamma-ray spectra. Nuclear               //
//   Instruments and Methods in Physics Research A 401 (1997) 113-         //
//   132.                                                                  //
//                                                                         //
//   [2]  M.Morhac et al.: Efficient one- and two-dimensional Gold         //
//   deconvolution and its application to gamma-ray spectra                //
//   decomposition. Nuclear Instruments and Methods in Physics             //
//   Research A 401 (1997) 385-408.                                        //
//                                                                         //
//   [3]  M.Morhac et al.: Identification of peaks in multidimensional     //
//   coincidence gamma-ray spectra. Nuclear Instruments and Methods in     //
//   Research Physics A  443(2000), 108-125.                               //
//                                                                         //
//   These NIM papers are also available as doc or ps files from:          //
//____________________________________________________________________________

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>


double       fResolution;     //resolution of the neighboring peaks
int         fgAverageWindow; //Average window of searched peaks
int         fgIterations;    //Maximum number of decon iterations (default=3)

   enum {
       kBackOrder2 =0,
       kBackOrder4 =1,
       kBackOrder6 =2,
       kBackOrder8 =3,
       kBackIncreasingWindow =0,
       kBackDecreasingWindow =1,
       kBackSmoothing3 =3,
       kBackSmoothing5 =5,
       kBackSmoothing7 =7,
       kBackSmoothing9 =9,
       kBackSmoothing11 =11,
       kBackSmoothing13 =13,
       kBackSmoothing15 =15
   };

int SpectrumfgIterations    = 3;
int SpectrumfgAverageWindow = 3;

#define PEAK_WINDOW 1024


/////////////////////NEW FUNCTIONS  JANUARY 2006
SEXP R_SpectrumBackground(SEXP R_spectrum,
                                          SEXP R_numberIterations,
                                          SEXP R_direction, SEXP R_filterOrder,
                                          SEXP R_smoothing,SEXP R_smoothWindow,
                                          SEXP R_compton)
{
  double * spectrum=REAL(R_spectrum);
  int numberIterations=INTEGER(R_numberIterations)[0];
  int ssize=LENGTH(R_spectrum);
  int direction=INTEGER(R_direction)[0];
  int filterOrder=INTEGER(R_filterOrder)[0];
  int smoothing=INTEGER(R_smoothing)[0];
  int smoothWindow=INTEGER(R_smoothWindow)[0];
  int compton=INTEGER(R_compton)[0];
  SEXP f;
/////////////////////////////////////////////////////////////////////////////
//        ONE-DIMENSIONAL BACKGROUND ESTIMATION FUNCTION - GENERAL FUNCTION
//
//        This function calculates background spectrum from source spectrum.
//        The result is placed in the vector pointed by spe1945ctrum pointer.
//
//        Function parameters:
//        spectrum-pointer to the vector of source spectrum
//        ssize-length of the spectrum vector
//        numberIterations-maximal width of clipping window,
//        direction- direction of change of clipping window
//               - possible values=kBackIncreasingWindow
//                                 kBackDecreasingWindow
//        filterOrder-order of clipping filter,
//                  -possible values=kBackOrder2
//                                   kBackOrder4
//                                   kBackOrder6
//                                   kBackOrder8
//        smoothing- logical variable whether the smoothing operation
//               in the estimation of background will be included
//             - possible values=FALSE
//                               TRUE
//        smoothWindow-width of smoothing window,
//                  -possible values=kBackSmoothing3
//                                   kBackSmoothing5
//                                   kBackSmoothing7
//                                   kBackSmoothing9
//                                   kBackSmoothing11
//                                   kBackSmoothing13
//                                   kBackSmoothing15
//         compton- logical variable whether the estimation of Compton edge
//                  will be included
//             - possible values=FALSE
//                               TRUE
//
///////////////////////////////////////////////////////////////////////////////
//

   int i, j, w, bw, b1, b2, priz;
   double a, b, c, d, e, yb1, yb2, ai, av, men, b4, c4, d4, e4, b6, c6, d6, e6, f6, g6, b8, c8, d8, e8, f8, g8, h8, i8;
   if (ssize <= 0)
      Rf_error ("Wrong Parameters");
   if (numberIterations < 1)
      Rf_error( "Width of Clipping Window Must Be Positive");
   if (ssize < 2 * numberIterations + 1)
      Rf_error( "Too Large Clipping Window");
   if (smoothing == TRUE && smoothWindow != kBackSmoothing3 && smoothWindow != kBackSmoothing5 && smoothWindow != kBackSmoothing7 && smoothWindow != kBackSmoothing9 && smoothWindow != kBackSmoothing11 && smoothWindow != kBackSmoothing13 && smoothWindow != kBackSmoothing15)
      Rf_error( "Incorrect width of smoothing window");
   double working_space[2 * ssize];
   for (i = 0; i < ssize; i++){
      working_space[i] = spectrum[i];
      working_space[i + ssize] = spectrum[i];
   }
   bw=(smoothWindow-1)/2;
   if (direction == kBackIncreasingWindow)
      i = 1;
   else if(direction == kBackDecreasingWindow)
      i = numberIterations;
   if (filterOrder == kBackOrder2) {
      do{
         for (j = i; j < ssize - i; j++) {
            if (smoothing == FALSE){
               a = working_space[ssize + j];
               b = (working_space[ssize + j - i] + working_space[ssize + j + i]) / 2.0;
               if (b < a)
                  a = b;
               working_space[j] = a;
            }

            else if (smoothing == TRUE){
               a = working_space[ssize + j];
               av = 0;
               men = 0;
               for (w = j - bw; w <= j + bw; w++){
                  if ( w >= 0 && w < ssize){
                     av += working_space[ssize + w];
                     men +=1;
                  }
               }
               av = av / men;
               b = 0;
               men = 0;
               for (w = j - i - bw; w <= j - i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     b += working_space[ssize + w];
                     men +=1;
                  }
               }
               b = b / men;
               c = 0;
               men = 0;
               for (w = j + i - bw; w <= j + i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     c += working_space[ssize + w];
                     men +=1;
                  }
               }
               c = c / men;
               b = (b + c) / 2;
               if (b < a)
                  av = b;
               working_space[j]=av;
            }
         }
         for (j = i; j < ssize - i; j++)
            working_space[ssize + j] = working_space[j];
         if (direction == kBackIncreasingWindow)
            i+=1;
         else if(direction == kBackDecreasingWindow)
            i-=1;
      } while ((direction == kBackIncreasingWindow && i <= numberIterations || direction == kBackDecreasingWindow && i >= 1));
   }

   else if (filterOrder == kBackOrder4) {
      do{
         for (j = i; j < ssize - i; j++) {
            if (smoothing == FALSE){
               a = working_space[ssize + j];
               b = (working_space[ssize + j - i] + working_space[ssize + j + i]) / 2.0;
               c = 0;
               ai = i / 2;
               c -= working_space[ssize + j - (int) (2 * ai)] / 6;
               c += 4 * working_space[ssize + j - (int) ai] / 6;
               c += 4 * working_space[ssize + j + (int) ai] / 6;
               c -= working_space[ssize + j + (int) (2 * ai)] / 6;
               if (b < c)
                  b = c;
               if (b < a)
                  a = b;
               working_space[j] = a;
            }

            else if (smoothing == TRUE){
               a = working_space[ssize + j];
               av = 0;
               men = 0;
               for (w = j - bw; w <= j + bw; w++){
                  if ( w >= 0 && w < ssize){
                     av += working_space[ssize + w];
                     men +=1;
                  }
               }
               av = av / men;
               b = 0;
               men = 0;
               for (w = j - i - bw; w <= j - i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     b += working_space[ssize + w];
                     men +=1;
                  }
               }
               b = b / men;
               c = 0;
               men = 0;
               for (w = j + i - bw; w <= j + i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     c += working_space[ssize + w];
                     men +=1;
                  }
               }
               c = c / men;
               b = (b + c) / 2;
               ai = i / 2;
               b4 = 0, men = 0;
               for (w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     b4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               b4 = b4 / men;
               c4 = 0, men = 0;
               for (w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     c4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               c4 = c4 / men;
               d4 = 0, men = 0;
               for (w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     d4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               d4 = d4 / men;
               e4 = 0, men = 0;
               for (w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     e4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               e4 = e4 / men;
               b4 = (-b4 + 4 * c4 + 4 * d4 - e4) / 6;
               if (b < b4)
                  b = b4;
               if (b < a)
                  av = b;
               working_space[j]=av;
            }
         }
         for (j = i; j < ssize - i; j++)
            working_space[ssize + j] = working_space[j];
         if (direction == kBackIncreasingWindow)
            i+=1;
         else if(direction == kBackDecreasingWindow)
            i-=1;
      }while(direction == kBackIncreasingWindow && i <= numberIterations || direction == kBackDecreasingWindow && i >= 1);
   }

   else if (filterOrder == kBackOrder6) {
      do{
         for (j = i; j < ssize - i; j++) {
            if (smoothing == FALSE){
               a = working_space[ssize + j];
               b = (working_space[ssize + j - i] + working_space[ssize + j + i]) / 2.0;
               c = 0;
               ai = i / 2;
               c -= working_space[ssize + j - (int) (2 * ai)] / 6;
               c += 4 * working_space[ssize + j - (int) ai] / 6;
               c += 4 * working_space[ssize + j + (int) ai] / 6;
               c -= working_space[ssize + j + (int) (2 * ai)] / 6;
               d = 0;
               ai = i / 3;
               d += working_space[ssize + j - (int) (3 * ai)] / 20;
               d -= 6 * working_space[ssize + j - (int) (2 * ai)] / 20;
               d += 15 * working_space[ssize + j - (int) ai] / 20;
               d += 15 * working_space[ssize + j + (int) ai] / 20;
               d -= 6 * working_space[ssize + j + (int) (2 * ai)] / 20;
               d += working_space[ssize + j + (int) (3 * ai)] / 20;
               if (b < d)
                  b = d;
               if (b < c)
                  b = c;
               if (b < a)
                  a = b;
               working_space[j] = a;
            }

            else if (smoothing == TRUE){
               a = working_space[ssize + j];
               av = 0;
               men = 0;
               for (w = j - bw; w <= j + bw; w++){
                  if ( w >= 0 && w < ssize){
                     av += working_space[ssize + w];
                     men +=1;
                  }
               }
               av = av / men;
               b = 0;
               men = 0;
               for (w = j - i - bw; w <= j - i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     b += working_space[ssize + w];
                     men +=1;
                  }
               }
               b = b / men;
               c = 0;
               men = 0;
               for (w = j + i - bw; w <= j + i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     c += working_space[ssize + w];
                     men +=1;
                  }
               }
               c = c / men;
               b = (b + c) / 2;
               ai = i / 2;
               b4 = 0, men = 0;
               for (w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     b4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               b4 = b4 / men;
               c4 = 0, men = 0;
               for (w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     c4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               c4 = c4 / men;
               d4 = 0, men = 0;
               for (w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     d4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               d4 = d4 / men;
               e4 = 0, men = 0;
               for (w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     e4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               e4 = e4 / men;
               b4 = (-b4 + 4 * c4 + 4 * d4 - e4) / 6;
               ai = i / 3;
               b6 = 0, men = 0;
               for (w = j - (int)(3 * ai) - bw; w <= j - (int)(3 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     b6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               b6 = b6 / men;
               c6 = 0, men = 0;
               for (w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     c6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               c6 = c6 / men;
               d6 = 0, men = 0;
               for (w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     d6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               d6 = d6 / men;
               e6 = 0, men = 0;
               for (w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     e6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               e6 = e6 / men;
               f6 = 0, men = 0;
               for (w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     f6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               f6 = f6 / men;
               g6 = 0, men = 0;
               for (w = j + (int)(3 * ai) - bw; w <= j + (int)(3 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     g6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               g6 = g6 / men;
               b6 = (b6 - 6 * c6 + 15 * d6 + 15 * e6 - 6 * f6 + g6) / 20;
               if (b < b6)
                  b = b6;
               if (b < b4)
                  b = b4;
               if (b < a)
                  av = b;
               working_space[j]=av;
            }
         }
         for (j = i; j < ssize - i; j++)
            working_space[ssize + j] = working_space[j];
         if (direction == kBackIncreasingWindow)
            i+=1;
         else if(direction == kBackDecreasingWindow)
            i-=1;
      }while(direction == kBackIncreasingWindow && i <= numberIterations || direction == kBackDecreasingWindow && i >= 1);
   }

   else if (filterOrder == kBackOrder8) {
      do{
         for (j = i; j < ssize - i; j++) {
            if (smoothing == FALSE){
               a = working_space[ssize + j];
               b = (working_space[ssize + j - i] + working_space[ssize + j + i]) / 2.0;
               c = 0;
               ai = i / 2;
               c -= working_space[ssize + j - (int) (2 * ai)] / 6;
               c += 4 * working_space[ssize + j - (int) ai] / 6;
               c += 4 * working_space[ssize + j + (int) ai] / 6;
               c -= working_space[ssize + j + (int) (2 * ai)] / 6;
               d = 0;
               ai = i / 3;
               d += working_space[ssize + j - (int) (3 * ai)] / 20;
               d -= 6 * working_space[ssize + j - (int) (2 * ai)] / 20;
               d += 15 * working_space[ssize + j - (int) ai] / 20;
               d += 15 * working_space[ssize + j + (int) ai] / 20;
               d -= 6 * working_space[ssize + j + (int) (2 * ai)] / 20;
               d += working_space[ssize + j + (int) (3 * ai)] / 20;
               e = 0;
               ai = i / 4;
               e -= working_space[ssize + j - (int) (4 * ai)] / 70;
               e += 8 * working_space[ssize + j - (int) (3 * ai)] / 70;
               e -= 28 * working_space[ssize + j - (int) (2 * ai)] / 70;
               e += 56 * working_space[ssize + j - (int) ai] / 70;
               e += 56 * working_space[ssize + j + (int) ai] / 70;
               e -= 28 * working_space[ssize + j + (int) (2 * ai)] / 70;
               e += 8 * working_space[ssize + j + (int) (3 * ai)] / 70;
               e -= working_space[ssize + j + (int) (4 * ai)] / 70;
               if (b < e)
                  b = e;
               if (b < d)
                  b = d;
               if (b < c)
                  b = c;
               if (b < a)
                  a = b;
               working_space[j] = a;
            }

            else if (smoothing == TRUE){
               a = working_space[ssize + j];
               av = 0;
               men = 0;
               for (w = j - bw; w <= j + bw; w++){
                  if ( w >= 0 && w < ssize){
                     av += working_space[ssize + w];
                     men +=1;
                  }
               }
               av = av / men;
               b = 0;
               men = 0;
               for (w = j - i - bw; w <= j - i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     b += working_space[ssize + w];
                     men +=1;
                  }
               }
               b = b / men;
               c = 0;
               men = 0;
               for (w = j + i - bw; w <= j + i + bw; w++){
                  if ( w >= 0 && w < ssize){
                     c += working_space[ssize + w];
                     men +=1;
                  }
               }
               c = c / men;
               b = (b + c) / 2;
               ai = i / 2;
               b4 = 0, men = 0;
               for (w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     b4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               b4 = b4 / men;
               c4 = 0, men = 0;
               for (w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     c4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               c4 = c4 / men;
               d4 = 0, men = 0;
               for (w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     d4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               d4 = d4 / men;
               e4 = 0, men = 0;
               for (w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     e4 += working_space[ssize + w];
                     men +=1;
                  }
               }
               e4 = e4 / men;
               b4 = (-b4 + 4 * c4 + 4 * d4 - e4) / 6;
               ai = i / 3;
               b6 = 0, men = 0;
               for (w = j - (int)(3 * ai) - bw; w <= j - (int)(3 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     b6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               b6 = b6 / men;
               c6 = 0, men = 0;
               for (w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     c6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               c6 = c6 / men;
               d6 = 0, men = 0;
               for (w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     d6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               d6 = d6 / men;
               e6 = 0, men = 0;
               for (w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     e6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               e6 = e6 / men;
               f6 = 0, men = 0;
               for (w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     f6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               f6 = f6 / men;
               g6 = 0, men = 0;
               for (w = j + (int)(3 * ai) - bw; w <= j + (int)(3 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     g6 += working_space[ssize + w];
                     men +=1;
                  }
               }
               g6 = g6 / men;
               b6 = (b6 - 6 * c6 + 15 * d6 + 15 * e6 - 6 * f6 + g6) / 20;
               ai = i / 4;
               b8 = 0, men = 0;
               for (w = j - (int)(4 * ai) - bw; w <= j - (int)(4 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     b8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               b8 = b8 / men;
               c8 = 0, men = 0;
               for (w = j - (int)(3 * ai) - bw; w <= j - (int)(3 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     c8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               c8 = c8 / men;
               d8 = 0, men = 0;
               for (w = j - (int)(2 * ai) - bw; w <= j - (int)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     d8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               d8 = d8 / men;
               e8 = 0, men = 0;
               for (w = j - (int)ai - bw; w <= j - (int)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     e8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               e8 = e8 / men;
               f8 = 0, men = 0;
               for (w = j + (int)ai - bw; w <= j + (int)ai + bw; w++){
                  if (w >= 0 && w < ssize){
                     f8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               f8 = f8 / men;
               g8 = 0, men = 0;
               for (w = j + (int)(2 * ai) - bw; w <= j + (int)(2 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     g8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               g8 = g8 / men;
               h8 = 0, men = 0;
               for (w = j + (int)(3 * ai) - bw; w <= j + (int)(3 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     h8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               h8 = h8 / men;
               i8 = 0, men = 0;
               for (w = j + (int)(4 * ai) - bw; w <= j + (int)(4 * ai) + bw; w++){
                  if (w >= 0 && w < ssize){
                     i8 += working_space[ssize + w];
                     men +=1;
                  }
               }
               i8 = i8 / men;
               b8 = ( -b8 + 8 * c8 - 28 * d8 + 56 * e8 - 56 * f8 - 28 * g8 + 8 * h8 - i8)/70;
               if (b < b8)
                  b = b8;
               if (b < b6)
                  b = b6;
               if (b < b4)
                  b = b4;
               if (b < a)
                  av = b;
               working_space[j]=av;
            }
         }
         for (j = i; j < ssize - i; j++)
            working_space[ssize + j] = working_space[j];
         if (direction == kBackIncreasingWindow)
            i += 1;
         else if(direction == kBackDecreasingWindow)
            i -= 1;
      }while(direction == kBackIncreasingWindow && i <= numberIterations || direction == kBackDecreasingWindow && i >= 1);
   }

   if (compton == TRUE) {
      for (i = 0, b2 = 0; i < ssize; i++){
         b1 = b2;
         a = working_space[i], b = spectrum[i];
         j = i;
         if (abs(a - b) >= 1) {
            b1 = i - 1;
            if (b1 < 0)
               b1 = 0;
            yb1 = working_space[b1];
            for (b2 = b1 + 1, c = 0, priz = 0; priz == 0 && b2 < ssize; b2++){
               a = working_space[b2], b = spectrum[b2];
               c = c + b - yb1;
               if (abs(a - b) < 1) {
                  priz = 1;
                  yb2 = b;
               }
            }
            if (b2 == ssize)
               b2 -= 1;
            yb2 = working_space[b2];
            if (yb1 <= yb2){
               for (j = b1, c = 0; j <= b2; j++){
                  b = spectrum[j];
                  c = c + b - yb1;
               }
               if (c > 1){
                  c = (yb2 - yb1) / c;
                  for (j = b1, d = 0; j <= b2 && j < ssize; j++){
                     b = spectrum[j];
                     d = d + b - yb1;
                     a = c * d + yb1;
                     working_space[ssize + j] = a;
                  }
               }
            }

            else{
               for (j = b2, c = 0; j >= b1; j--){
                  b = spectrum[j];
                  c = c + b - yb2;
               }
               if (c > 1){
                  c = (yb1 - yb2) / c;
                  for (j = b2, d = 0;j >= b1 && j >= 0; j--){
                     b = spectrum[j];
                     d = d + b - yb2;
                     a = c * d + yb2;
                     working_space[ssize + j] = a;
                  }
               }
            }
            i=b2;
         }
      }
   }
   PROTECT(f = allocVector(REALSXP,ssize));
   for (j = 0; j < ssize; j++){
      REAL(f)[j] = working_space[ssize + j];
   }
   UNPROTECT(1);
   return(f);
}


SEXP R_SpectrumSmoothMarkov(SEXP R_source, SEXP R_averWindow)
{
  double * source=REAL(R_source);
  int ssize=LENGTH(R_source);
  int averWindow=INTEGER(R_averWindow)[0];
  SEXP f;
/////////////////////////////////////////////////////////////////////////////
//        ONE-DIMENSIONAL MARKOV SPECTRUM SMOOTHING FUNCTION
//
//        This function calculates smoothed spectrum from source spectrum
//        based on Markov chain method.
//        The result is placed in the array pointed by source pointer.
//
//        Function parameters:
//        source-pointer to the array of source spectrum
//        ssize-length of source array
//        averWindow-width of averaging smoothing window
//
/////////////////////////////////////////////////////////////////////////////
   int xmin = 0, xmax= ssize - 1, i, l;
   double a, b, maxch;
   double nom, nip, nim, sp, sm, area = 0;
   if(averWindow <= 0)
      Rf_error( "Averaging Window must be positive");
   double working_space[ssize];
   for(i = 0, maxch = 0; i < ssize; i++){
      working_space[i]=0;
      if(maxch < source[i])
         maxch = source[i];

      area += source[i];
   }
   if(maxch == 0)
      return 0 ;

   nom = 1;
   working_space[xmin] = 1;
   for(i = xmin; i < xmax; i++){
      nip = source[i] / maxch;
      nim = source[i + 1] / maxch;
      sp = 0,sm = 0;
      for(l = 1; l <= averWindow; l++){
         if((i + l) > xmax)
            a = source[xmax] / maxch;

         else
            a = source[i + l] / maxch;
         b = a - nip;
         if(a + nip <= 0)
            a = 1;

         else
            a = sqrt(a + nip);
         b = b / a;
         b = exp(b);
         sp = sp + b;
         if((i - l + 1) < xmin)
            a = source[xmin] / maxch;

         else
            a = source[i - l + 1] / maxch;
         b = a - nim;
         if(a + nim <= 0)
            a = 1;
         else
            a = sqrt(a + nim);
         b = b / a;
         b = exp(b);
         sm = sm + b;
      }
      a = sp / sm;
      a = working_space[i + 1] = working_space[i] * a;
      nom = nom + a;
   }
   for(i = xmin; i <= xmax; i++){
      working_space[i] = working_space[i] / nom;
   }
   PROTECT(f = allocVector(REALSXP,ssize));
   for (i = 0; i < ssize; i++){
      REAL(f)[i] = working_space[i]*area;
   }
   UNPROTECT(1);
   return(f);

}


SEXP R_SpectrumDeconvolution(SEXP R_source, SEXP R_response,
                                      SEXP R_numberIterations,
                                      SEXP R_numberRepetitions, SEXP R_boost )
{

  double *source = REAL(R_source);
  double *response=REAL(R_response);
  int ssize=LENGTH(R_source);
  int numberIterations=INTEGER(R_numberIterations)[0];
  int numberRepetitions=INTEGER(R_numberRepetitions)[0];
  double boost=REAL(R_boost)[0];
  SEXP f;
/////////////////////////////////////////////////////////////////////////////
//   ONE-DIMENSIONAL DECONVOLUTION FUNCTION                                //
//   This function calculates deconvolution from source spectrum           //
//   according to response spectrum using Gold algorithm                   //
//   The result is placed in the vector pointed by source pointer.         //
//                                                                         //
//   Function parameters:                                                  //
//   source:  pointer to the vector of source spectrum                     //
//   response:     pointer to the vector of response spectrum              //
//   ssize:    length of source and response spectra                       //
//   numberIterations, for details we refer to the reference given below   //
//   numberRepetitions, for repeated boosted deconvolution                 //
//   boost, boosting coefficient                                           //
//                                                                         //
//    M. Morhac, J. Kliman, V. Matousek, M. Veselsk?, I. Turzo.:           //
//    Efficient one- and two-dimensional Gold deconvolution and its        //
//    application to gamma-ray spectra decomposition.                      //
//    NIM, A401 (1997) 385-408.                                            //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//

   if (ssize <= 0)
      Rf_error( "Wrong Parameters");

   if (numberRepetitions <= 0)
      Rf_error( "Wrong Parameters ");

       //   working_space-pointer to the working vector
       //   (its size must be 4*ssize of source spectrum)
   double working_space[4 * ssize];
   int i, j, k, lindex, posit = 0, lh_gold = -1, l, repet;
   double lda, ldb, ldc, area=0, maximum=0;
//read response vector
   for (i = 0; i < ssize; i++) {
      lda = response[i];
      if (lda != 0)
         lh_gold = i + 1;
      working_space[i] = lda;
      area += lda;
      if (lda > maximum) {
         maximum = lda;
         posit = i;
      }
   }
   if (lh_gold == -1)
      Rf_error( "ZERO RESPONSE VECTOR" );

//read source vector
   for (i = 0; i < ssize; i++)
      working_space[2 * ssize + i] = source[i];

// create matrix at*a and vector at*y
   for (i = 0; i < ssize; i++){
      lda = 0;
      for (j = 0; j < ssize; j++){
         ldb = working_space[j];
         k = i + j;
         if (k < ssize){
            ldc = working_space[k];
            lda = lda + ldb * ldc;
         }
      }
      working_space[ssize + i] = lda;
      lda = 0;
      for (k = 0; k < ssize; k++){
         l = k - i;
         if (l >= 0){
            ldb = working_space[l];
            ldc = working_space[2 * ssize + k];
            lda = lda + ldb * ldc;
         }
      }
      working_space[3 * ssize + i]=lda;
   }

// move vector at*y
   for (i = 0; i < ssize; i++){
      working_space[2 * ssize + i] = working_space[3 * ssize + i];
   }

//initialization of resulting vector
   for (i = 0; i < ssize; i++)
      working_space[i] = 1;

       //**START OF ITERATIONS**
   for (repet = 0; repet < numberRepetitions; repet++) {
      if (repet != 0) {
         for (i = 0; i < ssize; i++)
            working_space[i] = pow(working_space[i], boost);
      }
      for (lindex = 0; lindex < numberIterations; lindex++) {
         for (i = 0; i < ssize; i++) {
            if (working_space[2 * ssize + i] > 0.000001
                 && working_space[i] > 0.000001) {
               lda = 0;
               for (j = 0; j < lh_gold; j++) {
                  ldb = working_space[j + ssize];
                  if (j != 0){
                     k = i + j;
                     ldc = 0;
                     if (k < ssize)
                        ldc = working_space[k];
                     k = i - j;
                     if (k >= 0)
                        ldc += working_space[k];
                  }

                  else
                     ldc = working_space[i];
                  lda = lda + ldb * ldc;
               }
               ldb = working_space[2 * ssize + i];
               if (lda != 0)
                  lda = ldb / lda;

               else
                  lda = 0;
               ldb = working_space[i];
               lda = lda * ldb;
               working_space[3 * ssize + i] = lda;
            }
         }
         for (i = 0; i < ssize; i++)
            working_space[i] = working_space[3 * ssize + i];
      }
   }

//shift and write back resulting spectrum
PROTECT(f = allocVector(REALSXP,ssize));
   for (i = 0; i < ssize; i++) {
      lda = working_space[i];
      j = i + posit;
      j = j % ssize;
      REAL(f)[j] = lda*area;
   }
UNPROTECT(1);
   return(f);
}

SEXP R_SpectrumDeconvolutionRL(SEXP R_source, SEXP R_response,
                                      SEXP R_numberIterations,
                                      SEXP R_numberRepetitions, SEXP R_boost )
{

  double *source = REAL(R_source);
  double *response=REAL(R_response);
  int ssize=LENGTH(R_source);
  int numberIterations=INTEGER(R_numberIterations)[0];
  int numberRepetitions=INTEGER(R_numberRepetitions)[0];
  double boost=REAL(R_boost)[0];
  SEXP f;
/////////////////////////////////////////////////////////////////////////////
//   ONE-DIMENSIONAL DECONVOLUTION FUNCTION                                //
//   This function calculates deconvolution from source spectrum           //
//   according to response spectrum using Richardson-Lucy algorithm        //
//   The result is placed in the vector pointed by source pointer.         //
//                                                                         //
//   Function parameters:                                                  //
//   source:  pointer to the vector of source spectrum                     //
//   response:     pointer to the vector of response spectrum              //
//   ssize:    length of source and response spectra                       //
//   numberIterations, for details we refer to the reference given above   //
//   numberRepetitions, for repeated boosted deconvolution                 //
//   boost, boosting coefficient                                           //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//

   if (ssize <= 0)
      Rf_error( "Wrong Parameters");

   if (numberRepetitions <= 0)
      Rf_error( "Wrong Parameters");

       //   working_space-pointer to the working vector
       //   (its size must be 4*ssize of source spectrum)
   double working_space[4 * ssize];
   int i, j, k, lindex, posit, lh_gold, repet, kmin, kmax;
   double lda, ldb, ldc, maximum;
   lh_gold = -1;
   posit = 0;
   maximum = 0;

//read response vector
   for (i = 0; i < ssize; i++) {
      lda = response[i];
      if (lda != 0)
         lh_gold = i + 1;
      working_space[ssize + i] = lda;
      if (lda > maximum) {
         maximum = lda;
         posit = i;
      }
   }
   if (lh_gold == -1)
      Rf_error( "ZERO RESPONSE VECTOR");

//read source vector
   for (i = 0; i < ssize; i++)
      working_space[2 * ssize + i] = source[i];

//initialization of resulting vector
   for (i = 0; i < ssize; i++){
      if (i <= ssize - lh_gold)
         working_space[i] = 1;

      else
         working_space[i] = 0;

   }
       //**START OF ITERATIONS**
   for (repet = 0; repet < numberRepetitions; repet++) {
      if (repet != 0) {
         for (i = 0; i < ssize; i++)
            working_space[i] = pow(working_space[i], boost);
      }
      for (lindex = 0; lindex < numberIterations; lindex++) {
         for (i = 0; i <= ssize - lh_gold; i++){
            lda = 0;
            if (working_space[i] > 0){//x[i]
               for (j = i; j < i + lh_gold; j++){
                  ldb = working_space[2 * ssize + j];//y[j]
                  if (j < ssize){
                     if (ldb > 0){//y[j]
                        kmax = j;
                        if (kmax > lh_gold - 1)
                           kmax = lh_gold - 1;
                        kmin = j + lh_gold - ssize;
                        if (kmin < 0)
                           kmin = 0;
                        ldc = 0;
                        for (k = kmax; k >= kmin; k--){
                           ldc += working_space[ssize + k] * working_space[j - k];//h[k]*x[j-k]
                        }
                        if (ldc > 0)
                           ldb = ldb / ldc;

                        else
                           ldb = 0;
                     }
                     ldb = ldb * working_space[ssize + j - i];//y[j]*h[j-i]/suma(h[j][k]x[k])
                  }
                  lda += ldb;
               }
               lda = lda * working_space[i];
            }
            working_space[3 * ssize + i] = lda;
         }
         for (i = 0; i < ssize; i++)
            working_space[i] = working_space[3 * ssize + i];
      }
   }

//shift and write back resulting spectrum
PROTECT(f = allocVector(REALSXP,ssize));
   for (i = 0; i < ssize; i++) {
      lda = working_space[i];
      j = i + posit;
      j = j % ssize;
      REAL(f)[j] = lda;
   }
   UNPROTECT(1);
   return(f);
}

const char *SpectrumUnfolding(double *source,
                                               const double **respMatrix,
                                               int ssizex, int ssizey,
                                               int numberIterations,
                                               int numberRepetitions, double boost)
{
/////////////////////////////////////////////////////////////////////////////
//        ONE-DIMENSIONAL UNFOLDING FUNCTION
//        This function unfolds source spectrum
//        according to response matrix columns.
//        The result is placed in the vector pointed by source pointer.
//
//        Function parameters:
//        source-pointer to the vector of source spectrum
//        respMatrix-pointer to the matrix of response spectra
//        ssizex-length of source spectrum and # of columns of response matrix
//        ssizey-length of destination spectrum and # of rows of
//              response matrix
//        numberIterations, for details we refer to manual
//        Note!!! ssizex must be >= ssizey
/////////////////////////////////////////////////////////////////////////////
   int i, j, k, lindex, lhx = 0, repet;
   double lda, ldb, ldc, area;
   if (ssizex <= 0 || ssizey <= 0)
      return "Wrong Parameters";
   if (ssizex < ssizey)
      return "Sizex must be greater than sizey)";
   if (numberIterations <= 0)
      return "Number of iterations must be positive";
   double working_space[ssizex * ssizey + 2 * ssizey * ssizey + 4 * ssizex];

/*read response matrix*/
   for (j = 0; j < ssizey && lhx != -1; j++) {
      area = 0;
      lhx = -1;
      for (i = 0; i < ssizex; i++) {
         lda = respMatrix[j][i];
         if (lda != 0) {
            lhx = i + 1;
         }
         working_space[j * ssizex + i] = lda;
         area = area + lda;
      }
      if (lhx != -1) {
         for (i = 0; i < ssizex; i++)
            working_space[j * ssizex + i] /= area;
      }
   }
   if (lhx == -1)
      return ("ZERO COLUMN IN RESPONSE MATRIX");

/*read source vector*/
   for (i = 0; i < ssizex; i++)
      working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex + i] =
          source[i];

/*create matrix at*a + at*y */
   for (i = 0; i < ssizey; i++) {
      for (j = 0; j < ssizey; j++) {
         lda = 0;
         for (k = 0; k < ssizex; k++) {
            ldb = working_space[ssizex * i + k];
            ldc = working_space[ssizex * j + k];
            lda = lda + ldb * ldc;
         }
         working_space[ssizex * ssizey + ssizey * i + j] = lda;
      }
      lda = 0;
      for (k = 0; k < ssizex; k++) {
         ldb = working_space[ssizex * i + k];
         ldc =
             working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex +
                           k];
         lda = lda + ldb * ldc;
      }
      working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i] =
          lda;
   }

/*move vector at*y*/
   for (i = 0; i < ssizey; i++)
      working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex + i] =
          working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i];

/*create matrix at*a*at*a + vector at*a*at*y */
   for (i = 0; i < ssizey; i++) {
      for (j = 0; j < ssizey; j++) {
         lda = 0;
         for (k = 0; k < ssizey; k++) {
            ldb = working_space[ssizex * ssizey + ssizey * i + k];
            ldc = working_space[ssizex * ssizey + ssizey * j + k];
            lda = lda + ldb * ldc;
         }
         working_space[ssizex * ssizey + ssizey * ssizey + ssizey * i + j] =
             lda;
      }
      lda = 0;
      for (k = 0; k < ssizey; k++) {
         ldb = working_space[ssizex * ssizey + ssizey * i + k];
         ldc =
             working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex +
                           k];
         lda = lda + ldb * ldc;
      }
      working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i] =
          lda;
   }

/*move at*a*at*y*/
   for (i = 0; i < ssizey; i++)
      working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex + i] =
          working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i];

/*initialization in resulting vector */
   for (i = 0; i < ssizey; i++)
      working_space[ssizex * ssizey + 2 * ssizey * ssizey + i] = 1;

        /***START OF ITERATIONS***/
   for (repet = 0; repet < numberRepetitions; repet++) {
      if (repet != 0) {
         for (i = 0; i < ssizey; i++)
            working_space[ssizex * ssizey + 2 * ssizey * ssizey + i] = pow(working_space[ssizex * ssizey + 2 * ssizey * ssizey + i], boost);
      }
      for (lindex = 0; lindex < numberIterations; lindex++) {
         for (i = 0; i < ssizey; i++) {
            lda = 0;
            for (j = 0; j < ssizey; j++) {
               ldb =
                   working_space[ssizex * ssizey + ssizey * ssizey + ssizey * i + j];
               ldc = working_space[ssizex * ssizey + 2 * ssizey * ssizey + j];
               lda = lda + ldb * ldc;
            }
            ldb =
                working_space[ssizex * ssizey + 2 * ssizey * ssizey + 2 * ssizex + i];
            if (lda != 0) {
               lda = ldb / lda;
            }

            else
               lda = 0;
            ldb = working_space[ssizex * ssizey + 2 * ssizey * ssizey + i];
            lda = lda * ldb;
            working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i] = lda;
         }
         for (i = 0; i < ssizey; i++)
            working_space[ssizex * ssizey + 2 * ssizey * ssizey + i] =
                working_space[ssizex * ssizey + 2 * ssizey * ssizey + 3 * ssizex + i];
      }
   }

/*write back resulting spectrum*/
   for (i = 0; i < ssizex; i++) {
      if (i < ssizey)
         source[i] = working_space[ssizex * ssizey + 2 * ssizey * ssizey + i];

      else
         source[i] = 0;
   }
   return 0;
}

SEXP R_SpectrumSearchHighRes(SEXP R_source,
                                     SEXP R_sigma, SEXP R_threshold,
                                     SEXP  R_backgroundRemove, SEXP R_deconIterations,
                                     SEXP  R_markov, SEXP  R_averWindow)
{
     double *source=REAL(R_source);
     int ssize=LENGTH(R_source);
     double sigma=REAL(R_sigma)[0];
     double threshold=REAL(R_threshold)[0];
     int backgroundRemove=INTEGER(R_backgroundRemove)[0];
     int deconIterations=INTEGER(R_deconIterations)[0];
     int markov=INTEGER(R_markov)[0];
     int averWindow=INTEGER(R_averWindow)[0];
     int fMaxPeaks=ssize;
     int fNPeaks;
     double fPositionX[fMaxPeaks];
     SEXP destVector,f,ans,ans_names;

/////////////////////////////////////////////////////////////////////////////
//        ONE-DIMENSIONAL HIGH-RESOLUTION PEAK SEARCH FUNCTION
//        This function searches for peaks in source spectrum
//      It is based on deconvolution method. First the background is
//      removed (if desired), then Markov spectrum is calculated
//      (if desired), then the response function is generated
//      according to given sigma and deconvolution is carried out.
//
//        Function parameters:
//        source-pointer to the vector of source spectrum
//        destVector-pointer to the vector of resulting deconvolved spectrum     */
//        ssize-length of source spectrum
//        sigma-sigma of searched peaks, for details we refer to manual
//        threshold-threshold value in % for selected peaks, peaks with
//                amplitude less than threshold*highest_peak/100
//                are ignored, see manual
//      backgroundRemove-logical variable, set if the removal of
//                background before deconvolution is desired
//      deconIterations-number of iterations in deconvolution operation
//      markov-logical variable, if it is true, first the source spectrum
//             is replaced by new spectrum calculated using Markov
//             chains method.
//        averWindow-averanging window of searched peaks, for details
//                  we refer to manual (applies only for Markov method)
//
/////////////////////////////////////////////////////////////////////////////
//
   int i, j, numberIterations = (int)(7 * sigma + 0.5);
   double a, b, c;
   int k, lindex, posit, imin, imax, jmin, jmax, lh_gold, priz;
   double lda, ldb, ldc, area, maximum, maximum_decon;
   int xmin, xmax, l, peak_index = 0, size_ext = ssize + 2 * numberIterations, shift = numberIterations, bw = 2, w;
   double maxch;
   double nom, nip, nim, sp, sm, plocha = 0;
   double m0low=0,m1low=0,m2low=0,l0low=0,l1low=0,detlow,av,men;
   if (sigma < 1) {
      Rf_error("SearchHighRes", "Invalid sigma, must be greater than or equal to 1");
      return 0;
   }

   if(threshold<=0 || threshold>=100){
      Rf_error("SearchHighRes", "Invalid threshold, must be positive and less than 100");
      return 0;
   }

   j = (int) (5.0 * sigma + 0.5);
   if (j >= PEAK_WINDOW / 2) {
      Rf_error("SearchHighRes", "Too large sigma");
      return 0;
   }

   if (markov == TRUE) {
      if (averWindow <= 0) {
         Rf_error("SearchHighRes", "Averanging window must be positive");
         return 0;
      }
   }

   if(backgroundRemove == TRUE){
      if(ssize < 2 * numberIterations + 1){
         Rf_error("SearchHighRes", "Too large clipping window");
         return 0;
      }
   }

   k = (int) (2 * sigma+0.5);
   if(k >= 2){
      for(i = 0;i < k;i++){
         a = i,b = source[i];
         m0low += 1,m1low += a,m2low += a * a,l0low += b,l1low += a * b;
      }
      detlow = m0low * m2low - m1low * m1low;
      if(detlow != 0)
         l1low = (-l0low * m1low + l1low * m0low) / detlow;

      else
         l1low = 0;
      if(l1low > 0)
         l1low=0;
   }

   else{
      l1low = 0;
   }

   i = (int)(7 * sigma + 0.5);
   i = 2 * i;
   double working_space[7 * (ssize + i)];
   for (j=0;j<7 * (ssize + i);j++) working_space[j] = 0;
   for(i = 0; i < size_ext; i++){
      if(i < shift){
         a = i - shift;
         working_space[i + size_ext] = source[0] + l1low * a;
         if(working_space[i + size_ext] < 0)
            working_space[i + size_ext]=0;
      }

      else if(i >= ssize + shift){
      	 a = i - (ssize - 1 + shift);
         working_space[i + size_ext] = source[ssize - 1];
         if(working_space[i + size_ext] < 0)
            working_space[i + size_ext]=0;
      }

      else
         working_space[i + size_ext] = source[i - shift];
   }

   if(backgroundRemove == TRUE){
      for(i = 1; i <= numberIterations; i++){
         for(j = i; j < size_ext - i; j++){
            if(markov == FALSE){
               a = working_space[size_ext + j];
               b = (working_space[size_ext + j - i] + working_space[size_ext + j + i]) / 2.0;
               if(b < a)
                  a = b;

               working_space[j]=a;
            }

            else{
               a = working_space[size_ext + j];
               av = 0;
               men = 0;
               for (w = j - bw; w <= j + bw; w++){
                  if ( w >= 0 && w < size_ext){
                     av += working_space[size_ext + w];
                     men +=1;
                  }
               }
               av = av / men;
               b = 0;
               men = 0;
               for (w = j - i - bw; w <= j - i + bw; w++){
                  if ( w >= 0 && w < size_ext){
                     b += working_space[size_ext + w];
                     men +=1;
                  }
               }
               b = b / men;
               c = 0;
               men = 0;
               for (w = j + i - bw; w <= j + i + bw; w++){
                  if ( w >= 0 && w < size_ext){
                     c += working_space[size_ext + w];
                     men +=1;
                  }
               }
               c = c / men;
               b = (b + c) / 2;
               if (b < a)
                  av = b;
               working_space[j]=av;
            }
         }
         for(j = i; j < size_ext - i; j++)
            working_space[size_ext + j] = working_space[j];
      }
      for(j = 0;j < size_ext; j++){
         if(j < shift){
         	  a = j - shift;
         	  b = source[0] + l1low * a;
         	  if(b < 0)
         	     b = 0;
            working_space[size_ext + j] = b - working_space[size_ext + j];
         }

         else if(j >= ssize + shift){
         	  a = j - (ssize - 1 + shift);
         	  b = source[ssize - 1];
         	  if(b < 0)
         	     b = 0;
            working_space[size_ext + j] = b - working_space[size_ext + j];
         }

         else{
            working_space[size_ext + j] = source[j - shift] - working_space[size_ext + j];
         }
      }
      for(j = 0;j < size_ext; j++){
      	if(working_space[size_ext + j] < 0)
      	   working_space[size_ext + j] = 0;
      }
   }

   for(i = 0; i < size_ext; i++){
      working_space[i + 6*size_ext] = working_space[i + size_ext];
   }

   if(markov == TRUE){
      for(j = 0; j < size_ext; j++)
         working_space[2 * size_ext + j] = working_space[size_ext + j];
      xmin = 0,xmax = size_ext - 1;
      for(i = 0, maxch = 0; i < size_ext; i++){
         working_space[i] = 0;
         if(maxch < working_space[2 * size_ext + i])
            maxch = working_space[2 * size_ext + i];
         plocha += working_space[2 * size_ext + i];
      }
      if(maxch == 0) {
         return 0;
      }

      nom = 1;
      working_space[xmin] = 1;
      for(i = xmin; i < xmax; i++){
         nip = working_space[2 * size_ext + i] / maxch;
         nim = working_space[2 * size_ext + i + 1] / maxch;
         sp = 0,sm = 0;
         for(l = 1; l <= averWindow; l++){
            if((i + l) > xmax)
               a = working_space[2 * size_ext + xmax] / maxch;

            else
               a = working_space[2 * size_ext + i + l] / maxch;

            b = a - nip;
            if(a + nip <= 0)
               a=1;

            else
               a = sqrt(a + nip);

            b = b / a;
            b = exp(b);
            sp = sp + b;
            if((i - l + 1) < xmin)
               a = working_space[2 * size_ext + xmin] / maxch;

            else
               a = working_space[2 * size_ext + i - l + 1] / maxch;

            b = a - nim;
            if(a + nim <= 0)
               a = 1;

            else
               a = sqrt(a + nim);

            b = b / a;
            b = exp(b);
            sm = sm + b;
         }
         a = sp / sm;
         a = working_space[i + 1] = working_space[i] * a;
         nom = nom + a;
      }
      for(i = xmin; i <= xmax; i++){
         working_space[i] = working_space[i] / nom;
      }
      for(j = 0; j < size_ext; j++)
         working_space[size_ext + j] = working_space[j] * plocha;
      for(j = 0; j < size_ext; j++){
         working_space[2 * size_ext + j] = working_space[size_ext + j];
      }
      if(backgroundRemove == TRUE){
         for(i = 1; i <= numberIterations; i++){
            for(j = i; j < size_ext - i; j++){
               a = working_space[size_ext + j];
               b = (working_space[size_ext + j - i] + working_space[size_ext + j + i]) / 2.0;
               if(b < a)
                  a = b;
               working_space[j] = a;
            }
            for(j = i; j < size_ext - i; j++)
               working_space[size_ext + j] = working_space[j];
         }
         for(j = 0; j < size_ext; j++){
            working_space[size_ext + j] = working_space[2 * size_ext + j] - working_space[size_ext + j];
         }
      }
   }
//deconvolution starts
   area = 0;
   lh_gold = -1;
   posit = 0;
   maximum = 0;
//generate response vector
   for(i = 0; i < size_ext; i++){
      lda = (double)i - 3 * sigma;
      lda = lda * lda / (2 * sigma * sigma);
      j = (int)(1000 * exp(-lda));
      lda = j;
      if(lda != 0)
         lh_gold = i + 1;

      working_space[i] = lda;
      area = area + lda;
      if(lda > maximum){
         maximum = lda;
         posit = i;
      }
   }
//read source vector
   for(i = 0; i < size_ext; i++)
      working_space[2 * size_ext + i] = abs(working_space[size_ext + i]);
//create matrix at*a(vector b)
   i = lh_gold - 1;
   if(i > size_ext)
      i = size_ext;

   imin = -i,imax = i;
   for(i = imin; i <= imax; i++){
      lda = 0;
      jmin = 0;
      if(i < 0)
         jmin = -i;
      jmax = lh_gold - 1 - i;
      if(jmax > (lh_gold - 1))
         jmax = lh_gold - 1;

      for(j = jmin;j <= jmax; j++){
         ldb = working_space[j];
         ldc = working_space[i + j];
         lda = lda + ldb * ldc;
      }
      working_space[size_ext + i - imin] = lda;
   }
//create vector p
   i = lh_gold - 1;
   imin = -i,imax = size_ext + i - 1;
   for(i = imin; i <= imax; i++){
      lda = 0;
      for(j = 0; j <= (lh_gold - 1); j++){
         ldb = working_space[j];
         k = i + j;
         if(k >= 0 && k < size_ext){
            ldc = working_space[2 * size_ext + k];
            lda = lda + ldb * ldc;
         }

      }
      working_space[4 * size_ext + i - imin] = lda;
   }
//move vector p
   for(i = imin; i <= imax; i++)
      working_space[2 * size_ext + i - imin] = working_space[4 * size_ext + i - imin];
//initialization of resulting vector
   for(i = 0; i < size_ext; i++)
      working_space[i] = 1;
//START OF ITERATIONS
   for(lindex = 0; lindex < deconIterations; lindex++){
      for(i = 0; i < size_ext; i++){
         if(abs(working_space[2 * size_ext + i]) > 0.00001 && abs(working_space[i]) > 0.00001){
            lda=0;
            jmin = lh_gold - 1;
            if(jmin > i)
               jmin = i;

            jmin = -jmin;
            jmax = lh_gold - 1;
            if(jmax > (size_ext - 1 - i))
               jmax=size_ext-1-i;

            for(j = jmin; j <= jmax; j++){
               ldb = working_space[j + lh_gold - 1 + size_ext];
               ldc = working_space[i + j];
               lda = lda + ldb * ldc;
            }
            ldb = working_space[2 * size_ext + i];
            if(lda != 0)
               lda = ldb / lda;

            else
               lda = 0;

            ldb = working_space[i];
            lda = lda * ldb;
            working_space[3 * size_ext + i] = lda;
         }
      }
      for(i = 0; i < size_ext; i++){
         working_space[i] = working_space[3 * size_ext + i];
      }
   }
//shift resulting spectrum
   for(i=0;i<size_ext;i++){
      lda = working_space[i];
      j = i + posit;
      j = j % size_ext;
      working_space[size_ext + j] = lda;
   }
//write back resulting spectrum
   maximum = 0, maximum_decon = 0;
   j = lh_gold - 1;
   for(i = 0; i < size_ext - j; i++){
      if(i >= shift && i < ssize + shift){
         working_space[i] = area * working_space[size_ext + i + j];
         if(maximum_decon < working_space[i])
            maximum_decon = working_space[i];
         if(maximum < working_space[6 * size_ext + i])
            maximum = working_space[6 * size_ext + i];
      }

      else
         working_space[i] = 0;
   }
   lda=1;
   if(lda>threshold)
      lda=threshold;
   lda=lda/100;

//searching for peaks in deconvolved spectrum
   for(i = 1; i < size_ext - 1; i++){
      if(working_space[i] > working_space[i - 1] && working_space[i] > working_space[i + 1]){
         if(i >= shift && i < ssize + shift){
            if(working_space[i] > lda*maximum_decon && working_space[6 * size_ext + i] > threshold * maximum / 100.0){
               for(j = i - 1, a = 0, b = 0; j <= i + 1; j++){
                  a += (double)(j - shift) * working_space[j];
                  b += working_space[j];
               }
               a = a / b;
               if(a < 0)
                  a = 0;

               if(a >= ssize)
                  a = ssize - 1;
               if(peak_index == 0){
                  fPositionX[0] = a;
                  peak_index = 1;
               }

               else{
                  for(j = 0, priz = 0; j < peak_index && priz == 0; j++){
                     if(working_space[6 * size_ext + shift + (int)a] > working_space[6 * size_ext + shift + (int)fPositionX[j]])
                        priz = 1;
                  }
                  if(priz == 0){
                     if(j < fMaxPeaks){
                        fPositionX[j] = a;
                     }
                  }

                  else{
                     for(k = peak_index; k >= j; k--){
                        if(k < fMaxPeaks){
                           fPositionX[k] = fPositionX[k - 1];
                        }
                     }
                     fPositionX[j - 1] = a;
                  }
                  if(peak_index < fMaxPeaks)
                     peak_index += 1;
               }
            }
         }
      }
   }

   fNPeaks = peak_index;
   PROTECT(destVector = allocVector(REALSXP,ssize));
   for (i = 0; i < ssize; i++){
      REAL(destVector)[i] = working_space[shift+i];
   }
   UNPROTECT(1);
   PROTECT(f = allocVector(INTSXP,fNPeaks));
   for (i = 0; i < fNPeaks; i++){
     /*to account for 1-based vectros in R*/
      INTEGER(f)[i] = (int)fPositionX[i]+1;
   }
   UNPROTECT(1);
   PROTECT(ans = allocVector(VECSXP,2));
   PROTECT(ans_names = allocVector(VECSXP,2));
   SET_VECTOR_ELT(ans_names,1,Rf_mkString("y"));
   SET_VECTOR_ELT(ans_names,0,Rf_mkString("pos"));
   SET_VECTOR_ELT(ans,1,destVector);
   SET_VECTOR_ELT(ans,0,f);
   setAttrib(ans, R_NamesSymbol, ans_names);
   UNPROTECT(2);
   if(peak_index == fMaxPeaks)
      Rf_warning("SearchHighRes", "Peak buffer full");
   return(ans);
}

