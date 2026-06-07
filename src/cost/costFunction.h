//
// Filename     : costFunction.cc
// Description  : Classes defining cost functions
// Author(s)    : Pontus Melke (pontus@thep.lu.se)
// Created      : August 2006
// Revision     : $Id: costFunction.h 656 2016-01-21 14:33:38Z andre $
//

#ifndef COSTFUNCTION_H
#define COSTFUNCTION_H

#include "baseCostFunction.h"

///
/// @brief Cost function using a mean-square-error function
///
/// @details MeanSquare is defined as::
///
///@f[ \frac{1}{T}\frac{1}{N}\frac{1}{M}\sum_{t_1}^{T}\sum_{i=1}^{N}\sum_{j=1}^{M}(x_{ij}(t) - \tilde{x}_{ij}(t))^2, @f]
///
///where @f$T@f$ is the number of time points stored in the vector
///costTemplateTime_,@f$N@f$ is the number of compartments, @f$M@f$ is the
///number of species used in the cost calculation (specified in costList_),
///@f$x_{ij}(t)@f$ is the value of species @f$j@f$ in compartment @f$i@f$
///obtained from simulation and @f$\tilde{x}_{ij}(t)@f$ is the experimental
///value of that same quantity.
///
class MeanSquare : public BaseCostFunction {
 public:
	double getCost(double simValue, double dataValue, double scaleFactor);
};

///
/// @brief Cost function using a normalized mean-square-error function
///
/// MeanSquareRelative is defined as::
///
///@f[ \frac{1}{T}\frac{1}{N}\frac{1}{M}\sum_{t_1}^{T}\sum_{i=1}^{N}\sum_{j=1}^{M}\left(\frac{x_{ij}(t) - \tilde{x}_{ij}(t)}{\tilde{x}_{ij}(t)}\right)^2, @f]
///
///where @f$T@f$ is the number of time points,@f$N@f$ is the number of
///compartments, @f$M@f$ is the number of species used in the cost calculation
///(specified in costList_), @f$x_{ij}(t)@f$ is the value of species @f$j@f$
///in compartment @f$i@f$ obtained from simulation and @f$\tilde{x}_{ij}(t)@f$
///is the experimental value of that same quantity.
///
class MeanSquareRelative : public BaseCostFunction {
 public:
	double getCost(double simValue, double dataValue, double scaleFactor);
};

///
/// @brief Cost function using a mean-square-error function normalized
/// for different spatial regions and concentrations
///
/// NormalizedMeanSquare is defined as::
///
///@f[ \frac{1}{T}\frac{1}{M}\sum_{t_1}^{T}\sum_{j=1}^{M}\sum_{g=1}^{N_g}\frac{1}{N}\sum_{i_g=1}^{N}\left(\frac{x_{ij}(t) - \tilde{x}_{ij}(t)}{\tilde{x}_{ij}(t)}\right)^2, @f]
///
/// where @f$T@f$ is the number of time points,@f$N@f$ is the number of
/// compartments of a specific type (e.g. a clv region with 1 should be
/// normalized equally as the non clv3 region of many more cells), @f$M@f$
/// is the number of species used in the cost calculation
/// (specified in costList_), @f$x_{ij}(t)@f$ is the value of species @f$j@f$
/// in compartment @f$i@f$ obtained from simulation and @f$\tilde{x}_{ij}(t)@f$
/// is the experimental value of that same quantity.
///
class NormalizedMeanSquare : public BaseCostFunction {
 public:
  double getCost(double simValue, double dataValue, double scaleFactor);
  double addCost(const std::vector< std::vector<double> > &y, const double t); 
};

///
/// @brief Cost function using a max-square-error function useful for min-max optimizing,
///
///@f[ \max( \left( x_{ij}(t) - \tilde{x}_{ij}(t)\right)^2 )_{i \in M,j \in N, t \in T}, @f]
///
/// where @f$T@f$ is the set of different time points, @f$N@f$ is the set of different compartments,
/// @f$M@f$ is the set of species used in the cost.
/// @f$x_{ij}(t)@f$ is the value of species @f$j@f$ in compartment @f$i@f$ obtained from simulation
///  and @f$\tilde{x}_{ij}(t)@f$ is the experimental value of that same quantity.
///
class MaxSquare: public BaseCostFunction {
 public:
  double getCost(double simValue, double dataValue);
  double addCost(const std::vector< std::vector<double> > &y, const double t); 
};

///
/// @brief Cost function using a mean-square-error function normalized
/// for different spatial regions and concentrations. Works generically
///
/// DomainWeigthedMeanSquare is defined as::
///
///@f[ \frac{1}{T}\frac{1}{M}\frac{1}{N}\sum_{t}^{T}\sum_{m=1}^{M}\sum_{dm=1}^{D_m} \frac{N}{D_m N_{domain\_m}}\sum_{n=1}^{N_{domain\_m}}\left(x_{m,dm,n}(t) - \tilde{x}_{m,dm,n}(t)\right)^2, @f]
///
/// where @f$T@f$ is the number of time points,@f$N_{domain\_m}@f$ is the number of
/// compartments of a domain. @f$M@f$
/// is the number of species used in the cost calculation,
///  @f$x_{m, dm, n}(t)@f$ is the value of species @f$m@f$
/// in compartment @f$n@f$ obtained from simulation and @f$\tilde{x}_{m,dm,n}(t)@f$
/// is the experimental value of that same quantity.
///
///	The basic idea behind this is to increase the weight of compartments in a domain
/// so that they do not get run over by the non-domain compartment in the optimization.
/// If, for example, we are looking for an expression domain of only 1 compartment in
/// a template with 100 compartments where 99 of them are supposed to have zero value 
/// the cost-landscape traversed by the optimizer would in the original MeanSquare case 
/// be such that the optimizer would most likely find parameters such that all the 
/// compartments (the one in our sought after domain as well) get zero-valued.
/// This is avoided by saying that it is equally important that this single compartment
/// in the domain has the right concentration as it is that all the other compartments 
/// have zero value.
///
/// Domain identity of the compartments is decided by their value in the costtemplate file.
/// If the costtemplate specifies 0 for a compartment then this is not considered part of
/// the domain. All the compartments for which the value is non-zero is considered to be
/// part of the domain.
///


class DomainWeightedMeanSquare : public BaseCostFunction {
 public:
  double getCost(double simValue, double dataValue);
  double addCost(const std::vector< std::vector<double> > &y, const double t); 
};


///
/// @brief Cost function using a manhattan distance measure
///
/// @details The Manhattan cost is defined as::
///
///@f[ \frac{1}{T}\frac{1}{N}\frac{1}{M}\sum_{t_1}^{T}\sum_{i=1}^{N}\sum_{j=1}^{M} abs(x_{ij}(t) - \tilde{x}_{ij}(t)), @f]
///
/// where @f$T@f$ is the number of time points stored in the vector
/// costTemplateTime_,@f$N@f$ is the number of compartments, @f$M@f$ is the
/// number of species used in the cost calculation (specified in costList_),
/// @f$x_{ij}(t)@f$ is the value of species @f$j@f$ in compartment @f$i@f$
/// obtained from simulation and @f$\tilde{x}_{ij}(t)@f$ is the experimental
/// value of that same quantity.
///
class Manhattan : public BaseCostFunction {
 public:
  double getCost(double simValue, double dataValue, double scaleFactor);
};


///
/// @brief Cost function using cosine similarity measure
///
/// @details The cosine a between two vectors @f$ y, \tilde{y} @f$ are calculated according to
///
///@f[ cos(a) = \frac{y(t) \cdot \tilde{y}(t)}{|y_i| |\tilde{y}_i|}  @f]
///
/// The cosine similary is calculated for each pair of templates (simulation compared to data), and the average is taken for all templates.
/// This error should be low if the shape of the concentrations profiles of two templates are similar. The total cosine error is calculated as
///
///@f[ E_tot = \frac{1}{T}\sum_{t_1}^{T} \frac{y(t) \cdot \tilde{y}(t)}{|y| |\tilde{y}|}  @f]
///
/// where T is the number of time points. The simulation values of all species concentrations in all cells are stored 
/// in the vector @f$y@f$ while the experimental result is stored in the vector @f$\tilde{y}_i@f$.
/// @f$|y|@f$ indicates taking the norm of vector y, and @f$ x \cdot y@f$ is the dot product between vectors x and y.
///
class Cosine : public BaseCostFunction {
 public:
  double getCost(std::vector< double > simValue, std::vector< double > dataValue);
  double addCost(const std::vector< std::vector<double> > &y, const double t); 
};

/// @brief Cost function using a cell-cell cosine similarity measure
///
/// @details The cosine a between two vectors @f$ y, \tilde{y} @f$ are calculated according to
///
///@f[ cos(a) = \frac{y(t) \cdot \tilde{y}(t)}{|y_i| |\tilde{y}_i|}  @f]
///
/// The cosine similary is calculated for each pair of cells (simulation compared to data), and the average is taken for all included cells.
/// This error should be low if the shape of the concentrations profiles of two cells are similar. The total cosine error is calculated as
///
///@f[ E_tot = \frac{1}{T}\frac{1}{N}\sum_{t_1}^{T}\sum_{i=1}^{N} \frac{y_i(t) \cdot \tilde{y}_i(t)}{|y_i| |\tilde{y}_i|}  @f]
///
/// where N is the number of cells and T is the number of time points. The values of the species concentrations of  cell i from the
/// simulation result are stored as a vector@f$y_i@f$ while the experimental result for cell i is stored in the vector @f$\tilde{y}_i@f$.
/// @f$|y|@f$ indicates taking the norm of vector y, and @f$ x \cdot y@f$ is the dot product between vectors x and y.
///
class CosineCell : public BaseCostFunction {
 public:
  double getCost(std::vector< double > simValue, std::vector< double > dataValue);
  double addCost(const std::vector< std::vector<double> > &y, const double t); 
};



///
/// @brief Cost function using a species-species cosine similarity measure
///
/// @details The cosine a between two vectors @f$ y, \tilde{y} @f$ are calculated according to
///
///@f[ cos(a) = \frac{y(t) \cdot \tilde{y}(t)}{|y_i| |\tilde{y}_i|}  @f]
///
/// The cosine similary is calculated for each species  (simulation compared to data), and the average is taken for all
/// included cost species. This error should be low if the shape of the spatial pattern of a species is similar between
/// simulation/data. The full total cosine error between are calculated as
///
///@f[ E_tot = \frac{1}{T}\frac{1}{M}\sum_{t_1}^{T}\sum_{i=1}^{N} \frac{y_i(t) \cdot \tilde{y}_i(t)}{|y_i| |\tilde{y}_i|}  @f]
///
/// where M is the number of cells and T is the number of time points. The values of species i from the simulation result
/// are stored seen as a vector@f$y_i@f$ while the experimental result for species i is stored the vector @f$\tilde{y}_i@f$.
/// @f$|y|@f$ indicates taking the norm of vector y, and @f$ x \cdot y@f$ is the dot product between vectors x and y.
///
class CosineSpatial : public BaseCostFunction {
 public:
  double getCost(std::vector< double > simValue, std::vector< double > dataValue);
  double addCost(const std::vector< std::vector<double> > &y, const double t); 
};


///
/// @brief Composite meanSquare and cosineSpatial cost function
///
/// @details The cost is for each species calculated as
///
///@f[ E_{species} = w \cdot E_s + (1 - w) \cdot E_c  @f]
///
/// where E_s is the MeanSquare cost and E_c is the CosineSpatial cost (see documentation)
/// The weight w can only be set in the source code. This cost is then summed and averaged over all species/templates/time points
/// to give the total cost.
///
class  mscs: public BaseCostFunction {
 public:
  double getCost(std::vector< double > simValue, std::vector< double > dataValue);
  double addCost(const std::vector< std::vector<double> > &y, const double t); 
};


///
/// @brief Cost function based on relative entropy or Kullback-Leibler
///
/// @details Kullback-Leibler cost is defined as:
///
///@f[ \frac{1}{T}\frac{1}{N}\frac{1}{M}\sum_{t_1}^{T}\sum_{i=1}^{N}\sum_{j=1}^{M}( ln(\frac{\tilde{x}_{ij}(t)}{x_{ij}(t)}) \cdot \tilde{x}_{ij}(t))^2, @f]
///
/// where @f$T@f$ is the number of time points stored in the vector
/// costTemplateTime_,@f$N@f$ is the number of compartments, @f$M@f$ is the
/// number of species used in the cost calculation (specified in costList_),
/// @f$x_{ij}(t)@f$ is the value of species @f$j@f$ in compartment @f$i@f$
/// obtained from simulation and @f$\tilde{x}_{ij}(t)@f$ is the experimental
/// value of that same quantity. This error treats the simulation/experimental 
/// data sets as probability distributions and calculates the relative entropy
/// between them. Note that the simulation/experimental data sets are both
/// normalized such that the area-under-the-curve is one.
///
class KullbackLeibler : public BaseCostFunction {
 public:
	double getCost(double simValue, double dataValue);
	double addCost(const std::vector< std::vector<double> > &y, const double t); 
};




#endif
