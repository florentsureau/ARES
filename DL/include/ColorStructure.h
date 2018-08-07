/******************************************************************************
 **                   Copyright (C) 2016 by CEA
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author: Florent Sureau
 **
 **    Date:  06-07/2016
 **
 **    File:  ColorStructure.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Introduces class for Color Data
 **    -----------
 ** - 07/2017:  RGB and HSV color Spaces with conversion
 **
 ******************************************************************************/

#ifndef ColorStructure_H
#define ColorStructure_H

#include "DataStructure.h"

/**
 * @enum ColorSpace
 * @brief choice among color space for data
*/
enum class ColorSpace {RGB,HSV,LCH,LAB,HSV2,HSL};//FCS: only first two implemented

/**
 * @class ColorStruc
 * @brief Class describing the structure of color data considered as the product of
   different manifolds. Assumption is that a patch is either multichannel contiguous in a
   patch, i.e. for a RGB sample, we would have in one patch: all channel R, all channel G,
   all channel B and then these patches are concatenated across patches ; or the data is
   fully multichannel contiguous: we have all patches for a channel contiguous
   @see DataStruc::m_contiguousChannels flag.
*/
class ColorStruc : public DataStruc {
    private:
        ColorSpace m_cspace;
    public:
        ISAP::PatchExtractor3D* m_extractor;//!< Extractor to get patches
        ExtractorType m_extract_mode;//!< Type of patch extraction
        /**
          * Constructors.
        */
        ColorStruc(ISAP::PatchExtractor3D* _p,ExtractorType _m,size_t npatch=0,
                ColorSpace _cspace=ColorSpace::RGB):DataStruc(_p,_m,npatch){
            setColorSpace(_cspace);
        }

        ColorStruc(dblarray* inputData, size_t PatchSize2D,
                                        bool isContiguous,ColorSpace _cspace);

        /**
          * Destructor.
        */
        virtual ~ColorStruc(){ }

        /**
          * Signature of class.
        */
        virtual std::string getSignature() const ;

        //Getters



        //Setters

        /**
          * Set ColorSpace
           * @param[in] _c: colorspace for the data.
        */
        void setColorSpace(ColorSpace _c) { m_cspace=_c;}

        //Processing
        /**
          * Ensure RGB  ColorSpace Data
        **/
        void setToRGB();

        /**
          * Ensure HSV  ColorSpace Data
        **/
        void setToHSV();

        /**
          * Ensure LCH  ColorSpace Data
        **/
        void setToLCH();


    private:
        /**
          * RGB to HSV ColorSpace Data conversion
        **/
        void RGBtoHSV();
        /**
          * HSV to RGB ColorSpace Data conversion
        **/
        void HSVtoRGB();

        /**
          * RGB to HSV ColorSpace value conversion
           * @param[in] R: ptr to Red value.
           * @param[in] G: ptr to Green value.
           * @param[in] B: ptr to Blue value.
           * @param[out] H: ptr to Hue value.
           * @param[out] S: ptr to Saturation value.
           * @param[out] V: ptr to Brightness Value.
           @note this is a safe routine, i.e. H, S and V could be the same as
           R, G, B
        */
        inline void RGBtoHSV(double* R,double* G, double* B,
                            double* H, double* S, double* V){
            double maxRGB=MAX(MAX(*R,*G),*B);
            double minRGB=MIN(MIN(*R,*G),*B);
            double Chr=maxRGB-minRGB;
            //Convention in case all R,G,B are equal
            if(Chr!=0){
                if(maxRGB==*R) {
                    double diff=(*G - *B)/Chr;
                    *H=60 *(diff-6*floor(diff/6.0));
                }
                else if (maxRGB==*G) *H=60 * (*B - *R)/Chr +120;
                else *H=60 * ((*R - *G)/Chr)+240;
                *S=Chr/maxRGB;
            } else {
                *H=0;
                *S=0;
            }
            *V=maxRGB;
        }

        /**
          * HSV to RGB ColorSpace value conversion
           * @param[in] H: ptr to Hue value.
           * @param[in] S: ptr to Saturation value.
           * @param[in] V: ptr to Brightness Value.
           * @param[out] R: ptr to Red value.
           * @param[out] G: ptr to Green value.
           * @param[out] B: ptr to Blue value.
           @note this is a safe routine, i.e. H, S and V could be the same as
           R, G, B
        */
        inline void HSVtoRGB(double* H,double* S, double* V,
                             double* R, double* G, double* B){
            double Chr=*V * *S;
            double nH=fmod(*H,360)/60;
            double X=Chr*(1 - fabs(fmod(nH,2) -1));
            double off=*V - Chr;
            //Convention in case all R,G,B are equal
            if(Chr!=0){
                if(nH<1){ *R=Chr; *G=X; *B=0;}
                else if (nH<2){ *R=X; *G=Chr; *B=0;}
                else if (nH<3){ *R=0; *G=Chr; *B=X;}
                else if (nH<4){ *R=0; *G=X; *B=Chr;}
                else if (nH<5){ *R=X; *G=0; *B=Chr;}
                else{ *R=Chr; *G=0;*B=X;}
            } else{ *R=0; *G=0; *B=0;}
            *R+=off;
            *G+=off;
            *B+=off;
        }
};



#endif //ColorStructure_H
