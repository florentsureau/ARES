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
 **    Date:  01-02/2017
 **
 **    File:  GenForwardBackward.h
 **
 *******************************************************************************
 **
 **    DESCRIPTION: Headers for Generalized Forward Backward
 **    -----------
 ** - 01-02/2017: class for Generalized Forward Backward
 **
 ******************************************************************************/

#ifndef GenForwardBackward_H
#define GenForwardBackward_H

#include "SparseCoding.h"
#include "S1Decomposition.h"
#include <vector>
#include <list>
#include <algorithm>
#include <functional>


class GenFBSparse_Settings {
protected:
    double m_epsilon;/**< Precision used for comparisons.*/
    size_t m_nit_gFB;/**< Number of genFB iterations.*/
    unsigned int m_maxSparse;/**< Maximal sparsity index.*/
    double  m_maxL2Res;/**< Maximal L2 norm if sparsity < Maximal sparsity.*/
    double m_lambda;/**< Overall lagrangian multiplier.*/
    double m_relax;/**< Relaxation parameter for intermediate variables.*/
    double m_w_gFB;/**< Weight of constraint vs sparsity in [0,1].*/
    bool m_affHullProj;/**< Add affine constraint on the code*/
    bool m_convHullProj;/**< Add simplex constraint on the code*/
    bool m_posSparse;/**< Add positivity constraint on the code*/
    bool m_reweight;/**< Reweighting should be perform*/
    bool m_hardthresh;/**< Hard Tresholding instead of soft*/
    bool m_decay;/**< Weight Decay*/
    bool m_singleProjGrad;/**< Use sparse projection onto simplex in one step*/
    dblarray* m_weightDecay;/**< contains hyperparameters for weights*/
    dblarray* m_reweightParams;/**< contains hyperparameters for reweighting*/
    bool m_verbose;/**< Verbosity flag.*/
    bool m_leftprecond;/**< Use Left preconditioning flag.*/
    Manifold m_manifold;/**< Manifold flag.*/
    //Note: posSparse + affine should be same as convex
public:
    GenFBSparse_Settings (){ init();}
    GenFBSparse_Settings (const GenFBSparse_Settings& other){
        init();
        *this=other;
    }
    ~GenFBSparse_Settings() {
        if(m_reweightParams!=nullptr) delete m_reweightParams;
        if(m_weightDecay!=nullptr) delete m_weightDecay;
     }
    /**
      * Set Verbosity.
       * @param[in] _verb flag for Verbosity.
    */
    void setVerbose(bool _verb){m_verbose=_verb;}
    /**
      * Set Manifold we are working on.
       * @param[in] _man Manifold flag.
    */
    void setManifold(Manifold _man){m_manifold = _man;}
    /**
      * Set Hard thresholding instead of soft.
       * @param[in] _hard Hard thresholding flag.
    */
    void setHardThr(bool _hard){m_hardthresh = _hard;}
    /**
      * Set Maximal sparsity tolerated.
       * @param[in] _maxsp maximal l0 norm of the code.
    */
    void setMaxSparsity(unsigned int _maxsp){m_maxSparse = _maxsp;}
    /**
      * Set maximal residual energy tolerated (if sparsity<maximal sparsity).
       * @param[in] _maxl2 maximal l2 norm of residual.
    */
    void setMaxL2Res(unsigned int _maxl2){m_maxL2Res = _maxl2;}
    /**
      * Use left preconditioning to get code.
       * @param[in] _leftp Left preconditioning flag.
    */
    void setLeftPrecond(bool _leftp){m_leftprecond = _leftp;}
    /**
      * Set epsilon used for float equality, comparisons.
       * @param[in] _eps precision value.
    */
    void setEpsilon(double _eps){m_epsilon = _eps;}
    /**
      * Set number of iterations for generalized Forward Backward.
       * @param[in] _nit number of iterations of FB algorithm.
    */
    void setNitgFB(const int _nit){m_nit_gFB =_nit;}
    /**
      * Set generalized FB sparse hyperparameter.
       * @param[in] _lambda hyperparameter weighting the l1 norm.
    */
    void setLambda(double _lambda){m_lambda = _lambda;}
    /**
      * Set single Projection coupled with gradient descent.
       * @param[in] _singproj projectino of gradient descent.
    */
    void setSingleProjGrad(double _singproj){m_singleProjGrad = _singproj;}
    /**
      * Set relaxation parameter for Generalized Forward Backward.
       * @param[in] _relax precision value.
    */
    void setRelax(double _relax){
        if((_relax>0)&&(_relax<=1)) m_relax = _relax;
        else {
            std::cout<<"Relax parameter should be >0 and <=1)"<<std::endl;
            MyExit::exit(EXIT_FAILURE);
        }
    }
    /**
      * Set weight related to sparsity intermediate variable.
    */
    void setReweight(bool _reweight){m_reweight = _reweight;}
    /**
     * Set reweight hyperparameter.
      * @param[in] _reweight ptr to array containing the decreasing sequence of
        hyperparameters used (continuation).
     */
    void setReweightHyper(dblarray* _reweight){
        if(m_reweightParams!=nullptr) delete m_reweightParams;
        m_reweightParams=nullptr;
        if(_reweight!=nullptr){
            m_reweightParams=new dblarray(*_reweight);
            m_reweight=true;
        }
    }
    /**
     * Set weight hyperparameters.
      * @param[in] _weight ptr to array containing the decreasing sequence of
        hyperparameters used (continuation).
     */
    void setWeightHyper(dblarray* _weight){
        if(m_weightDecay!=nullptr){
            std::cout<<"DEL STUFF"<<std::endl;
            delete m_weightDecay;
        }
        m_weightDecay =nullptr;
        if(_weight!=nullptr) {
            std::cout<<"ASSIGN STUFF"<<std::endl;
            m_weightDecay = new dblarray(*_weight);
            m_decay=true;
        }
    }
    /**
      * Set weight related to sparsity intermediate variable.
    */
    void setSparseWeight(double _weight){m_w_gFB=_weight;}
    /**
      * Set Affine hull projection.
       * @param[in] _aff flag for affine hull projection.
    */
    void setAffHullProj(bool _aff){m_affHullProj =_aff;}
    /**
      * Set Convex hull projection.
       * @param[in] _conv flag for convex hull projection.
    */
    void setConvHullProj(bool _conv){m_convHullProj= _conv;}
    /**
     * Set Positivity projection when applying sparsity constraint.
     * @param[in] _conv flag for convex hull projection.
     */
    void setPosSparseProj(bool _pos){m_posSparse= _pos;}

    /**
      * Get verbosity parameter.
       * @return verbosity.
    */
    bool getVerbose() const{return m_verbose;}

    /**
      * Return flag for left preconditioning to get code.
       * @return flag for left preconditioning.
    */
    bool getLeftPrecond() const{return m_leftprecond;}
    /**
      * Get Manifold.
       * @return Manifold we are working on.
    */
    Manifold getManifold() const{return m_manifold;}
    /**
      * Get Maximal sparsity tolerated.
       * @return maximal sparsity.
    */
    unsigned int getMaxSparsity() const{return m_maxSparse;}
    /**
      * Get Maximal residual error energy if maximal sparsity not reached.
       * @return minimal sparsity.
    */
    double getMaxL2Res() const{return m_maxL2Res;}

    /**
      * Get precision for floating point comparisons.
       * @return precision used.
    */
    double getEpsilon() const{return m_epsilon;}
    /**
      * Get number of iterations in gFB.
       * @return number of iterations.
    */
    size_t getNitgFB() const{return m_nit_gFB;}
    /**
      * Get threshold for sparsity constraint.
       * @return soft threshold value.
    */
    double getLambda() const{return m_lambda;}
    /**
      * Get relaxation parameter for gFB.
       * @return relaxation parameter.
    */
    double getRelax() const{return m_relax;}
    /**
      * Get weight related to sparsity intermediate variable.
       * @return weight for sparse term (in convex combination).
    */
    double getSparseWeight() const{return m_w_gFB;}
    /**
      * Get affine hull projection flag gFB.
       * @return affine hull projection flag.
    */
    bool getAffHullProj() const{return m_affHullProj;}
    /**
      * Get convex hull projection flag gFB.
       * @return convex hull projection flag.
    */
    bool getConvHullProj() const{return m_convHullProj;}

    /**
     * Get convex hull projection flag gFB.
      * @return convex hull projection flag.
     */
    bool getPosSparseProj() const{return m_posSparse;}

    /**
     * Get reweight flag.
      * @return true if reweighting should be performed.
     */
    bool getReweight() const{return m_reweight;}

    /**
     * Get Hard Thresholding flag.
      * @return true if hard thresholding instead of soft should be performed.
     */
    bool getHardThresh() const{return m_hardthresh;}

    /**
     * Get Decay flag.
      * @return true if decaying lambda should be performed.
     */
    bool getDecay() const {return m_decay;}
    /**
     * Get Number of steps for lambda decay flag.
      * @return number of lambdas.
     */
    unsigned int getNDecay() const {
        if(m_weightDecay ==nullptr) return (unsigned int) 0;
        else return (unsigned int) m_weightDecay->nx();
    }

    /**
     * Get reweight hyperparameter.
      * @return ptr to array containing the decreasing sequence of
        hyperparameters used (continuation).
     */
    dblarray* getReweightHyper() const{return m_reweightParams;}

    /**
     * Get weight hyperparameters.
      * @return ptr to array containing the decreasing sequence of
        hyperparameters used (continuation).
     */
    dblarray* getWeightHyper() const{return m_weightDecay;}

    /**
     * Get reweight hyperparameter for an iteration.
      * @param[in] it iteration to get the hyperparameter
      * @return hyperparameter weight or -1 if it fails.
     */
    double getWeightHyper(unsigned int it) const{
        if(it<getNDecay()) return (*m_weightDecay)(it);
        else return -1;
    }
    /**
     * Get reweight hyperparameter for an iteration.
      * @param[in] it iteration to get the hyperparameter
      * @return hyperparameter weight or -1 if it fails.
     */
    double getReweightHyper(unsigned int it) const{
        if(it<getNReweight()) return (*m_reweightParams)(it);
        else return -1;
    }
    /**
     * Get number of reweighting desired.
      * @return number of successive reweighting desired.
     */
    unsigned int getNReweight() const{
        if(m_reweightParams ==nullptr) return (unsigned int) 0;
        else return (unsigned int) m_reweightParams->nx();
    }
   /**
      * Get flag for single projection after gradient descent .
       * @return .
    */
    bool getSingleProjGrad() const{return m_singleProjGrad;}

    /**
      * Initialize gFB parameters.
    */
    void init();

    /**
      * Display gFB parameters.
    */
    void print() const;

    GenFBSparse_Settings& operator =(const GenFBSparse_Settings& _sett){
        m_verbose=_sett.getVerbose();
        m_epsilon=_sett.getEpsilon();
        m_nit_gFB=_sett.getNitgFB();
        m_lambda=_sett.getLambda();
        m_relax=_sett.getRelax();
        m_w_gFB=_sett.getSparseWeight();
        m_affHullProj=_sett.getAffHullProj();
        m_convHullProj =_sett.getConvHullProj();
        m_posSparse =_sett.getPosSparseProj();
        m_reweight =_sett.getReweight();
        if(_sett.getReweightHyper() !=nullptr){
            if(m_reweightParams!=nullptr) delete m_reweightParams;
            m_reweightParams=new dblarray(*(_sett.getReweightHyper()));
        } else m_reweightParams=nullptr;
        m_leftprecond= _sett.getLeftPrecond();
        m_manifold= _sett.getManifold();
        m_maxL2Res =_sett.getMaxL2Res();
        m_maxSparse=_sett.getMaxSparsity();
        m_hardthresh=_sett.getHardThresh();
        m_decay=_sett.getDecay();
        m_singleProjGrad =_sett.getSingleProjGrad();
        if(_sett.getWeightHyper() !=nullptr){
            if(m_weightDecay!=nullptr) delete m_weightDecay;
            m_weightDecay =new dblarray(*(_sett.getWeightHyper()));
        } else m_weightDecay =nullptr;
        return *this;
    }
};


/**
 * @class constGenFBSparse
 * @brief Class implementing generalized forward backward for sparse coding,
   with affine or convex constraints on the code.
*/
class constGenFBSparse : public SparseCoding{
    protected:
    GenFBSparse_Settings m_gFBsettings;/**< Parameters for gFB.*/
    //Stopping criteria
    double m_deltaStop;/**<Stopping criterion - each update l2 norm on intermediate
    variables is less than this value.*/
    dblarray *m_deltaInter;/**<Update l2 norm for each intermediate variable.*/
    //Cost variables
    bool m_trackCost;/**<flag to indicate we want to track the cost.*/
    dblarray* m_saveCost;
    dblarray* m_saveAffCost;
    std::string m_costFileName;/**<Fits FileName where to save cost function.*/
    std::string m_affFileName;/**<Fits FileName where to save affine distance.*/
    //Iterations and (Re)weighting variables
    unsigned int m_it;/**<Iteration number (in subproblem if reweighting).*/
    unsigned int m_itrew;/**< current Reweighting Iteration number.*/
    unsigned int m_itdecay;/**< current weighting Iteration number.*/
    unsigned int m_nitrew;/**< Number of reweighting iterations.*/
    unsigned int m_ndecay;/**< Number of decaying iterations.*/
    dblarray* m_sparseWeights;/**< Weights for sparse penalty.*/
    bool m_weightSparse;/**< Flag set if using (re)weighting.*/
    bool m_isWeightLocal;/**< Flag set if weights locally allocated.*/
    //Independent processing variables
    bool m_indepComputing;/**< Flag set if each input processed independently.*/
    size_t m_nSamp;/**< Either 1 (fully independent) or m_ninput.*/

public:
    //constructors/destructors
    /**
      * Standart constructor.
    */
    constGenFBSparse(GenFBSparse_Settings& _sett): SparseCoding(){
        m_gFBsettings.init();
        setParams(_sett);
        init();
        m_gFBsettings.setManifold(Manifold::Rn);
    }

    /**
      * Initialize inner variables.
    */
    void init();

    /**
      * Standart destructor.
    */
    virtual ~constGenFBSparse(){
        if(m_deltaInter!=nullptr) delete m_deltaInter;
        if((m_isWeightLocal)&&(m_sparseWeights!=nullptr))
                                                        delete m_sparseWeights;
        if(m_saveCost!=nullptr) delete m_saveCost;
        if(m_saveAffCost!=nullptr) delete m_saveAffCost;
    }

    //setters
    /**
      * Set parameters of gFB.
       * @param[in] _sett parameters.
    */
    virtual void setParams(GenFBSparse_Settings _sett){m_gFBsettings=_sett;}

    /**
      * Set parallel computing parameter.
       * @param[in] _indep true if independent computing for each sample.
    */
    void setIndepComputing(bool _indep){m_indepComputing = _indep;}

    /**
      * Set (diagonal) weights for sparsity term.
       * @param[in] _sparseWeights (diagonal) weights used for sparsity.
    */
    void setWeights(dblarray* _sparseWeights){
        if((m_isWeightLocal)&&(m_sparseWeights!=nullptr))
                                                        delete m_sparseWeights;
        m_sparseWeights = _sparseWeights;
        m_isWeightLocal=false;
    }

    /**
      * Set dictionary, overriding base class implementation.
       * @param[in] dictionary dictionary to set.
    */
    virtual void setDictionary(const dblarray &dictionary);

    /**
      * Set input metric, overriding base class implementation.
       * @param[in] metric input metric to set.
    */
    virtual void setInputMetric(const dblarray &metric);
    /**
      * Set code metric, overriding base class implementation.
       * @param[in] metric code metric to set.
    */
    virtual void setCodeMetric(const dblarray &metric);
    /**
      * Set matrix decomposition, overriding base class implementation.
       * @param[in] mdecomp decomposition to use.
    */
    virtual void setMatDecomp(MatrixDecomposition &mdecomp);
    /**
      * Set matrix decomposition, overriding base class implementation.
    */
    virtual void setMatDecomp();

    //getters
    /**
      * Return gFB parameters.
       * @brief call GenFBSparse_Settings::print.
    */
    void getParam(){m_gFBsettings.print();}

    /**
      * Return manifold we are working on.
    */
    virtual Manifold getManifold() const { return m_gFBsettings.getManifold();}
    /**
      * Set parallel computing parameter.
       * @param[in] _indep true if independent computing for each sample.
    */
    bool getIndepComputing() const{return m_indepComputing;}

    /**
      * Return current iteration number.
       * @return current iteration number.
    */
    unsigned int getCurIt(){return m_it+m_itrew* m_gFBsettings.getNitgFB();}
    /**
      * Interface for signature.
       * @return signature of constGenFBSparse
    */
    virtual std::string getSignature(){
        std::string Algo;
        if(m_gFBsettings.getConvHullProj())
            Algo=std::string("genForwardBackwardConvexHull");
        else if(m_gFBsettings.getAffHullProj())
            Algo=std::string("genForwardBackwardAffineHull");
        else Algo=std::string("(gen)ForwardBackward");
        if(m_gFBsettings.getPosSparseProj()) Algo.append("PosSparse");
        else Algo.append("Sparse");
        return Algo;
    }

    virtual std::string getName() { return std::string("gFB"); }

    //i/os
    /**
      * Save the convergence parameters of intermediate vars in fits file.
       * @param[in] costName name of fits file where to save cost function.
       * @param[in] AffName name of fits file where to save distance from affine
            hull.
    */
    void setCost(std::string& costName,std::string& AffName){
        m_trackCost=true;
        m_costFileName = costName;
        m_affFileName = AffName;
    }

    /**
      * Save a fits file for debugging.
       * @param[in] fitsName prefix of the filename.
       * @param[in] kit iteration number.
       * @param[in] data data to save.
    */
    void saveDebugFile(std::string fitsName, unsigned int kit, dblarray* data){
        std::stringstream strstr;
        strstr.str("");
        strstr<<kit;
        fitsName.append(strstr.str());
        strstr.clear();
        fitsName.append(".fits");
        fits_write_dblarr(fitsName.c_str(), *data);
    }

    //processing functions

    /**
      * Apply square root of weights to the data and dictionary.
       * @param[in] Input data to weight.
       * @param[in] Weights metric used for WLS.
       * @return input weighted by square root of weights
       * @note the dictionary is also weighted accordingly
    */

    dblarray* applySqrtWeights(dblarray*& Input, dblarray& Weights) {
        std::cout<<"NOT IMPLEMENTED YET"<<std::endl;
        return new dblarray(*Input);
    }

    /**
      * Get closest atom from a point.
       * @param[in] InVec point to consider.
       * @param[in] Dico dictionary buffer (row major, atoms as columns).
       * @param[in] Npix number of pixels in input point and for atoms.
       * @param[in] NatomsDico number of atoms in dictionary.
       * @param[in] RefL2 distance in between atoms and dictionary.
    */
    virtual size_t getClosestAtom(double* InVec,double* Dico,size_t Npix
                                            ,size_t NatomsDico, double& RefL2);


    /**
     * Get an input sample from its index .
       * @param[in] InData input Data
       * @param[in] sampleNb index of the sample
       * @param[out] InSample ptr to buffer where the result is stored.
       * @param[in] Channel index of channel to get
    */
     virtual void getInputSample(DataStruc*& InData,size_t sampleNb,
                                                        double *InSample);
    /**
      * Basic tests to check  internal variables.
       * @brief exits if basic tests fail.
    */
    void checkConsistency();

    /**
      * Basic tests to check inputs and internal variables.
       * @param[in] inputVector input data
       * @brief runs internal tests, checks consistency with input, and exits if
         basic tests fail.
    */
    void checkInputs(dblarray& inputVector);

    /**
      * Initialize cost array.
       * @param[in] NRew iterations for reweighting.
       * @param[in] NDecay iterations for lambda weighting.
       * @return cost array initialized.
    */
    dblarray* initCost(size_t NRew,size_t NDecay);

    /**
      * Compute the l2 norm of the residual.
       * @param[in] residual ptr to array of residuals.
       * @return l2 norm of residual for each row.
    */
    dblarray* computel2Cost(dblarray* residual);

    /**
      * Compute the l2 norm of the residual for one vector.
       * @param[in] residual ptr to residual.
       * @param[in] Npix number of pixels in residual.
       * @return l2 norm of residual.
    */
   double computel2Cost(double* residual,size_t Npix);

    /**
      * Compute the l1 norm of the code.
       * @param[in] code ptr to array of codes.
       * @param[in] sthr soft threshold used in algorithm (combination of
         Lipschitz, weighting of the constraints and sparse parameter).
       * @return l1 norm of code for each row.
    */
    dblarray* computel1Cost(dblarray* code, const double sthr);

    /**
      * Compute the l1 norm of one code.
       * @param[in] code ptr to code.
       * @param[in] soft threshold used in algorithm (combination of Lipschitz,
         weighting of the constraints and sparse parameter).
       * @param[in] index index of the current row processed
       * @param[in] Ncoefs number of coefficients in code.
       * @param[in] lambda sparse parameter
       * @return l1 norm of code for each row.
    */
    double computel1Cost(double* Code,const double sthr,const size_t index,
                                    const size_t Ncoefs,const double lambda);

    /**
      * Compute the (l2) distance of the codes to the affine hull.
       * @param[in] code ptr to array of codes.
       * @return euclidean distance of the codes to the affine hull constraint.
    */
    dblarray* computeAffl2(dblarray* code);

    /**
      * Compute the (l2) distance of the code to the affine hull.
       * @param[in] code ptr to code.
       * @param[in] Ncoefs number of coefficients in code.
       * @return euclidean distance of the code to the affine hull constraint.
    */
    double computeAffl2(double* Code, const size_t Ncoefs);

    /**
      * Get Lipschitz constant for gradient, using power method.
       * @param[in] xsubi seed for initialization of the power method.
       * @param[in] nit_max maximum number of iterations for power method.
       * @param[in] Matrix Matrix used in Lipschitz constant computation.
    */
    double computeLipschitzCst(unsigned short xsubi[3],int nit_max,
                                                             dblarray* Matrix);

    /**
      * Implementation to get codes for input examples using gFB.
       * @param[in] inputVector ptr to input set (1 sample per row).
       * @return input set code .
    */
    virtual dblarray* computeCode(DataStruc* inputVector);

    /**
      * Interface to get the sparse approximation from input codes.
       * @param[in] Code ptr to input code.
       * @param[in] Input ptr to input structure.
       * @return ptr to array with sparse approximation.
    */
    virtual dblarray* computeApprox(dblarray* Code,
                                            DataStruc* inputVector=nullptr){
         dblarray* Approx=m_matdecomp->matMultTransp(*Code,*m_dictionary);
         //delete Code;
         return Approx;
      };

      /**
        * Interface to get the sparse approximation from input codes.
         * @param[in] Code ptr to input code.
         * @param[out] Approx ptr to output sparse approximation.
         * @param[in] Input ptr to input structure.
         * @return ptr to array with sparse approximation.
      */
    virtual void computeApprox(dblarray* Code, dblarray* Approx,
                                                DataStruc* inputVector=nullptr){
           m_matdecomp->matMultTransp(*Approx,*Code,*m_dictionary);
    };


    /**
      * Get verbosity parameter.
       * @return verbosity.
    */
    bool isVerbose() const{return m_gFBsettings.getVerbose();};

protected:

    /**
      * Check whether cost function need to be computed at that iteration.
       * @param[in] kit iteration number.
       * @return true if cost function needed at that iteration.
    */
    bool computeCost();

    void searchLambdaIndepCode(double* InVec, size_t Npix,double mu,
                        const size_t index, double* Code, dblarray* Dico,
                                                dblarray* Weights);


    /**
      * Implementation to get the code of an input example using gFB.
       * @param[in] inputVector ptr to input sample.
       * @param[in] Npix number of pixels in the current sample.
       * @param[in] mu relaxation parameter for gradient descent.
       * @oaram[in] index index of the current row processed.
       * @oaram[in] Dico dictionary used for sparse coding (i.e. can be a
       restriction of the owned dictionary).
       * @oaram[out] Weights array storing the weights used in sparse coding.
       * @param[out] Code ptr to buffer where code will be stored.
    */
    virtual void computeIndepCode(double* InVec,const size_t Npix,
                            const double mu, const size_t index, double* Code,
                               dblarray* Dico,dblarray* Weights, double lambda,
                               double* SaveBuffer, double* SaveAffBuffer);

    /**
      * Implementation to get all codes for input examples using gFB.
       * @param[in] inputVector ptr to input set (1 sample per row).
       * @param[in] mu relaxation parameter for gradient descent.
       * @oaram[in] Dico dictionary used for sparse coding (i.e. can be a
       restriction of the owned dictionary).
       * @oaram[out] Weights array storing the weights used in sparse coding.
       * @return input set code .
    */
    virtual dblarray* computeAllCodes(dblarray* inputVector,const double mu,
                                             dblarray* Dico,dblarray* Weights);


    /**
      * Implementation to compute projection on affine k-sparse codes
       * @param[in,out] outBuf ptr to buffer with the current intermediate
          variable with sparsity constrain.
       * @param[in] xBuf ptr to buffer storing the current full term.
       * @param[in] Ncoefs number of coefs in vector to update.
       * @return projection onto the intersection of affine constraint and
       k-sparse codes.
    */
    void computeProjGSHP(double* outBuf,double* xBuf, const size_t Ncoefs);

    /**
      * Implementation to compute projection on k-sparse codes and simplex
       * @param[in,out] outBuf ptr to buffer with the current intermediate
          variable with sparsity constrain.
       * @param[in] xBuf ptr to buffer storing the current full term.
       * @param[in] Ncoefs number of coefs in vector to update.
       * @return projection onto the intersection of simplex and k-sparse codes.
    */
    void computeProjGSSP(double* outBuf,double* xBuf, const size_t Ncoefs);

    /**
      * Soft Thresholding Operator for Generalized forward backward.
       * @param[in] upBuf ptr to buffer storing the constant term (2x-Grad(f(x)).
       * @param[in,out] outBuf ptr to buffer with the current intermediate
          variable with sparsity constrain.
       * @param[in] xBuf ptr to buffer storing the current full term.
       * @param[in] sthr threshold for soft thresholding.
       * @param[in,out] delta update l2 norm for sparse variable.
     * @param[in] Ncoefs number of coefs in vector to update.
    */
    void computeSparseUpdate(double* upBuf,double* outBuf, double* xBuf,
                                 double sthr,double* delta,const size_t Ncoefs);

    /**
      * Weighted Soft Thresholding Operator for Generalized forward backward.
       * @param[in] upBuf ptr to buffer storing the constant term (2x-Grad(f(x)).
       * @param[in,out] outBuf ptr to buffer with the current intermediate
          variable with sparsity constrain.
       * @param[in] xBuf ptr to buffer storing the current full term.
       * @param[in] wbuf array of threshold per coefficient.
       * @param[in,out] delta update l2 norm for sparse variable.
     * @param[in] Ncoefs number of coefs in vector to update.
    */
    void computeWSparseUpdate(double* upBuf,double* outBuf,double* xBuf,
                                double* wBuf,double* delta,const size_t Ncoefs);

    /**
     * Asymetric Soft Thresholding Operator for Generalized forward backward.
     * @param[in] upBuf ptr to buffer storing the constant term (2x-Grad(f(x)).
     * @param[in,out] outBuf ptr to buffer with the current intermediate
     variable with sparsity+positivity constraint.
     * @param[in] xBuf ptr to buffer storing the current full term.
     * @param[in] sthr threshold for soft thresholding.
     * @param[in,out] delta update l2 norm for sparse variable.
     * @param[in] Ncoefs number of coefs in vector to update.
     */
    void computePositiveSparseUpdate(double* upBuf,double* outBuf, double* xBuf,
                                 double sthr,double* delta,const size_t Ncoefs);

    /**
     * Asymetric Weighted Soft Thresholding Operator for gen. forward backward.
      * @param[in] upBuf ptr to buffer storing the constant term (2x-Grad(f(x)).
      * @param[in,out] outBuf ptr to buffer with the current intermediate
        variable with sparsity+positivity constraint.
      * @param[in] xBuf ptr to buffer storing the current full term.
      * @param[in] wbuf array of threshold per coefficient.
      * @param[in,out] delta update l2 norm for sparse variable.
      * @param[in] Ncoefs number of coefs in vector to update.
     */
    void computeWPosSparseUpdate(double* upBuf,double* outBuf,double* xBuf,
                              double* wBuf,double* delta, const size_t Ncoefs);

    /**
      * Projection onto affine hull of atoms for Generalized forward backward.
       * @param[in] upBuf ptr to buffer storing  constant term (2x-Grad(f(x)).
       * @param[in,out] outBuf ptr to buffer with the current intermediate
          variable for the code with affine constraint.
       * @param[in] xBuf ptr to buffer storing the current full term.
       * @param[in,out] delta update l2 norm for const variable.
       * @param[in] Ncoefs number of coefs in vector to update.
    */
    void computeAffineHullUpdate(double* upBuf,double* outBuf, double* xBuf,
                                            double* delta, const size_t Ncoefs);

    /**
      * Projection onto convex hull of atoms for Generalized forward backward.
       * @param[in] upBuf ptr to buffer storing  constant term (2x-Grad(f(x)).
       * @param[in,out] outBuf ptr to buffer with the current intermediate
          variable for the code with cobvex constraint.
       * @param[in] xBuf ptr to buffer storing the current full term.
       * @param[in,out] delta update l2 norm for const variable.
       * @param[in] Ncoefs number of coefs in vector to update.
    */
    void computeConvexHullUpdate(double* upBuf,double* outBuf,double* xBuf,
                                            double* delta, const size_t Ncoefs);

    /**
      * Compute the residual based on current codes and the dictionary.
       * @param[in] Code current codes.
       * @param[in] InVector input array to compute the residual.
       * @return ptr array to residual.
    */
    virtual dblarray* computeResidual(dblarray& Code,dblarray& InVector,
                                                                dblarray& Dico);

    /**
      * Compute the residual based on a current code and the dictionary.
       * @param[in] Code current code.
       * @param[in] InVec input vector to compute the residual.
       * @param[out] Res ptr array to residual.
    */
    virtual void computeResidual(double* Code,double* InVec,double*& Res,
                                                                dblarray& Dico);

    /**
      * Compute the update with the gradient term in gFB.
       * @param[in] Code current codes.
       * @param[in] Residual current residuals.
       * @param[in] mu relaxation term.
       * @return ptr to array containing the update with the gradient term.
    */
    virtual dblarray* computeUpGrad(dblarray& Code, dblarray& Residual,
                                                    double mu, dblarray& Dico);

    /**
      * Compute the update with the gradient in gFB for a single input.
       * @param[in] Code current code.
       * @param[in] Res current residual.
       * @param[in] mu relaxation term.
       * @param[out] Upcode ptr to gradient term.
    */
    virtual void computeUpGrad(double* Code,double* Res,double mu,
                                               double*& UpCode, dblarray& Dico);

    /**
      * Compute the residual based on a current code and the dictionary.
       * @param[in] z1Buf array to (pos) sparsity component
       * @param[in] z2Buf array to affine or convex projection componen
       * @param[out] CodeBuf array containing the output code.
       * @param[in] Ncoefs number of coefs in the code.
       * @return ptr array to residual.
    */
    virtual void UpdateCode(double* z1Buf,double* z2Buf,double* CodeBuf,
                                                            const size_t Ncoefs);

    /**
      * Update Weights for reweighting strategies.
       * @param[in] Code array containing the current codes.
       * @param[out] Weight array containing the output weights.
       * @param[in] hyper hyperparameter to use for reweighting.
       * @param[in] sthr soft thresholding value.
    */
    void updateWeights(dblarray& Code, dblarray& Weight, const double hyper,
                                                            const double sthr);

    /**
      * Update Weights for reweighting strategies.
       * @param[in] Code array containing the current code.
       * @param[out] Weight array containing the output weights.
       * @param[in] hyper hyperparameter to use for reweighting.
       * @param[in] sthr soft thresholding value.
       * @oaram[in] index index of the current row processed.
       * @param[in] Ncoefs number of coefs in the code.
    */
    void updateWeights(double* Code, double* Weight, const double hyper,
                                      const double sthr, const size_t Ncoefs);


};

/**
 * @class constGenFBSparseS1
 * @brief class for generalized forward backward for sparse coding in S1,
   with affine or convex constraints on the code.
*/
class constGenFBSparseS1 : public constGenFBSparse{
protected:
    S1Decomposition* m_s1decomp;/**< Ptr to S1 operations.*/
    bool m_maxCvxBall;/**< Flag to specify if atoms are required to
    belong to a convex ball around current point [i.e. riemannian distance in
    [-pi/2, pi/2] for each pixel*/
    bool m_recomputeLipschitz;/**< Flag if Lipschitz need be recomputed.*/
    std::vector<unsigned int> m_atomInMaxCvxBall;/**< List of atoms inside the
        strongly convex ball around current point in S1^N.*/
    std::vector<unsigned int> m_sampleNoAtom;
    /**< List of atoms with no atom in
    the strongly convex ball.*/
public:

    //constructors/destructors
    /**
      * Standart constructor.
    */
    constGenFBSparseS1(GenFBSparse_Settings& _sett,bool _cvxBall,
                                bool _recLip=true): constGenFBSparse(_sett){
        init();
        setMaxCvxBall(_cvxBall);
        setRecompLip(_recLip);
        m_gFBsettings.setManifold(Manifold::S1n);
    }

    /**
      * Initialize inner variables.
    */
    void init();

    /**
      * Interface for signature.
       * @return signature of constGenFBSparse
    */
    virtual std::string getSignature(){
        std::string Algo;
        if(m_gFBsettings.getConvHullProj())
            Algo=std::string("S1genForwardBackwardConvexHull");
        else if(m_gFBsettings.getAffHullProj())
            Algo=std::string("S1genForwardBackwardAffineHull");
        else Algo=std::string("S1(gen)ForwardBackward");
        if(m_gFBsettings.getPosSparseProj()) Algo.append("PosSparse");
        else Algo.append("Sparse");
        if(getMaxCvxBall()) Algo.append("MaxCvxBall");
        return Algo;
    }

    /**
      * Standart destructor.
    */
    virtual ~ constGenFBSparseS1(){
        if(!m_sampleNoAtom.empty()) m_sampleNoAtom.clear();
        if(! m_atomInMaxCvxBall.empty()) m_atomInMaxCvxBall.clear();
    }


    //setters
    /**
      * Set if atoms looked in strongly convex ball around samples.
       * @param[in] _cvxBall flag set to true if strongly convex ball used.
    */
    virtual void setMaxCvxBall(bool _cvxBall){m_maxCvxBall = _cvxBall;}
    /**
      * Set if Lipschitz constant needs to be computed per input sample.
       * @param[in] _recLip flag set to true if Lipschitz constant per input
        sample.
    */
    virtual void setRecompLip(bool _recLip){m_recomputeLipschitz = _recLip;}

    //getters
    /**
      * Get strongly convex ball flag.
       * @return strongly convex ball flag.
    */
    bool getMaxCvxBall() const{ return m_maxCvxBall;}

    //Processing functions

    /**
      * Get closest atom from a point.
       * @param[in] InVec point to consider.
       * @param[in] Dico dictionary buffer (row major, atoms as columns).
       * @param[in] Npix number of pixels in input point and for atoms.
       * @param[in] NatomsDico number of atoms in dictionary.
       * @param[in] RefL2 distance in between atoms and dictionary.
    */
    virtual size_t getClosestAtom(double* InVec,double*Dico,size_t Npix
                                            ,size_t NatomsDico, double& RefL2);

    /**
      * Get the list of atoms in the dictionary that belong to a strongly convex
        ball around current vector of pixels in S1^N.
       * @param[in] InBuffer ptr to current vector of pixels.
       * @return input set code .
    */

    void getAtomsInMaxCvxBall(double* InBuffer);

    /**
      * Implementation to get codes for input examples using gFB on S1^N data.
       * @param[in] inputVector ptr to input set (1 sample per row).
       * @return input set code .
    */
    virtual dblarray* computeCode(DataStruc* inputVector);

    /**
      * Compute S1 Log matrix containing coordinates of atoms on a local normal
        chart.
       * @param[in] input pointer to reference angles for the tangent plane
      * @return matrix containing the coordinates of the atoms in local normal
        chart with origin given by the angles in input
     */
    dblarray* constructLogMat(double* input);

    /**
      * Implementation to get approximation for input examples using gFB on
      S1 data.
       * @param[in] inputVector ptr to input set (1 sample per row).
       * @return input set sparse approximation .
    */
    virtual dblarray* computeApprox(dblarray* Code,
                                            DataStruc* inputVector=nullptr);


    /**
      * Implementation to get the sparse approximation using gFB on S1 data.
       * @param[in] Code ptr to input code.
       * @param[out] Approx ptr to output sparse approximation.
       * @param[in] Input ptr to input structure.
       * @return ptr to array with sparse approximation.
    */
    virtual void computeApprox(dblarray* Code, dblarray* Approx,
         DataStruc* inputVector=nullptr);

};

#endif //GenForwardBackward_H
