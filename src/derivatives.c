/*
    Copyright (C) 2015 Tomas Flouri, Diego Darriba

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll.h"

static int sumtable_tipinner(pll_partition_t * partition,
                             unsigned int parent_clv_index,
                             double * parent_clv_vector,
                             unsigned int child_clv_index,
                             double * child_clv_vector,
                             const unsigned int * parent_scaler,
                             const unsigned int * child_scaler,
                             const unsigned int * params_indices,
                             double *sumtable)
{
  int retval;
  unsigned int i;
  unsigned int tip_clv_index;
  unsigned int inner_clv_index;
  unsigned int sites = partition->sites;
  const unsigned int * scaler;

  double ** eigenvecs = (double **)malloc(partition->rate_cats *
                                          sizeof(double *));
  double ** inv_eigenvecs = (double **)malloc(partition->rate_cats *
                                              sizeof(double *));
  double ** freqs = (double **)malloc(partition->rate_cats *
                                      sizeof(double *));
  if (!eigenvecs || !inv_eigenvecs || !freqs)
  {
    if (eigenvecs) free(eigenvecs);
    if (inv_eigenvecs) free(inv_eigenvecs);
    if (freqs) free(freqs);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  /* ascertaiment bias correction */
  if (partition->asc_bias_alloc)
    sites += partition->states;

  for (i = 0; i < partition->rate_cats; ++i)
  {
    eigenvecs[i] = partition->eigenvecs[params_indices[i]];
    inv_eigenvecs[i] = partition->inv_eigenvecs[params_indices[i]];
    freqs[i] = partition->frequencies[params_indices[i]];
  }

  /* find which of the two child nodes is the tip */
  if (parent_clv_index < partition->tips)
  {
    tip_clv_index = parent_clv_index;
    inner_clv_index = child_clv_index;
    scaler = child_scaler;
  }
  else
  {
    tip_clv_index = child_clv_index;
    inner_clv_index = parent_clv_index;
    scaler = parent_scaler;
  }

  retval = pll_core_update_sumtable_ti(partition->states,
                                       sites,
                                       partition->rate_cats,
                                       (inner_clv_index == parent_clv_index) ? parent_clv_vector : child_clv_vector,
                                       partition->tipchars[tip_clv_index],
                                       scaler,
                                       eigenvecs,
                                       inv_eigenvecs,
                                       freqs,
                                       partition->tipmap,
                                       partition->maxstates,
                                       sumtable,
                                       partition->attributes);

  free(freqs);
  free(eigenvecs);
  free(inv_eigenvecs);

  return retval;
}

static int sumtable_innerinner(pll_partition_t * partition,
                                unsigned int parent_clv_index,
                                double * parent_clv_vector,
                                unsigned int child_clv_index,
                                double * child_clv_vector,
                                const unsigned int * parent_scaler,
                                const unsigned int * child_scaler,
                                const unsigned int * params_indices,
                                double *sumtable)
{
  int retval;
  unsigned int i;
  unsigned int sites = partition->sites;

  double ** eigenvecs = (double **)malloc(partition->rate_cats *
                                          sizeof(double *));
  double ** inv_eigenvecs = (double **)malloc(partition->rate_cats *
                                              sizeof(double *));
  double ** freqs = (double **)malloc(partition->rate_cats *
                                      sizeof(double *));
  if (!eigenvecs || !inv_eigenvecs || !freqs)
  {
    if (eigenvecs) free(eigenvecs);
    if (inv_eigenvecs) free(inv_eigenvecs);
    if (freqs) free(freqs);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  /* ascertaiment bias correction */
  if (partition->asc_bias_alloc)
    sites += (unsigned int) partition->asc_additional_sites;

  for (i = 0; i < partition->rate_cats; ++i)
  {
    eigenvecs[i] = partition->eigenvecs[params_indices[i]];
    inv_eigenvecs[i] = partition->inv_eigenvecs[params_indices[i]];
    freqs[i] = partition->frequencies[params_indices[i]];
  }

  retval = pll_core_update_sumtable_ii(partition->states,
                                       sites,
                                       partition->rate_cats,
                                       parent_clv_vector,
                                       child_clv_vector,
                                       parent_scaler,
                                       child_scaler,
                                       eigenvecs,
                                       inv_eigenvecs,
                                       freqs,
                                       sumtable,
                                       partition->attributes);

  free(freqs);
  free(eigenvecs);
  free(inv_eigenvecs);

  return retval;
}

static int sumtable_repeats(pll_partition_t * partition,
                                unsigned int parent_clv_index,
                                double * parent_clv_vector,
                                unsigned int child_clv_index,
                                double * child_clv_vector,
                                const unsigned int * parent_scaler,
                                const unsigned int * child_scaler,
                                const unsigned int * params_indices,
                                double *sumtable)
{
  int retval;
  unsigned int i;
  unsigned int sites = partition->sites;

  double ** eigenvecs = (double **)malloc(partition->rate_cats *
                                          sizeof(double *));
  double ** inv_eigenvecs = (double **)malloc(partition->rate_cats *
                                              sizeof(double *));
  double ** freqs = (double **)malloc(partition->rate_cats *
                                      sizeof(double *));
  if (!eigenvecs || !inv_eigenvecs || !freqs)
  {
    if (eigenvecs) free(eigenvecs);
    if (inv_eigenvecs) free(inv_eigenvecs);
    if (freqs) free(freqs);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  /* ascertaiment bias correction */
  if (partition->asc_bias_alloc)
    sites += (unsigned int) partition->asc_additional_sites;

  for (i = 0; i < partition->rate_cats; ++i)
  {
    eigenvecs[i] = partition->eigenvecs[params_indices[i]];
    inv_eigenvecs[i] = partition->inv_eigenvecs[params_indices[i]];
    freqs[i] = partition->frequencies[params_indices[i]];
  }

  const unsigned int * parent_site_id =
    pll_get_site_id(partition, parent_clv_index);
  const unsigned int * child_site_id = 
    pll_get_site_id(partition, child_clv_index);

  unsigned int parent_ids = pll_get_sites_number(partition, parent_clv_index);
  unsigned int child_ids = pll_get_sites_number(partition, child_clv_index);
  unsigned int inv = parent_ids > child_ids;
  retval =
    pll_core_update_sumtable_repeats(partition->states,
                        sites,
                        inv ? child_ids : parent_ids,
                        partition->rate_cats,
                        inv ? child_clv_vector : parent_clv_vector,
                        !inv ? child_clv_vector : parent_clv_vector,
                        inv ? child_scaler : parent_scaler,
                        !inv ? child_scaler : parent_scaler,
                        eigenvecs,
                        inv_eigenvecs,
                        freqs,
                        sumtable,
                        inv ? child_site_id : parent_site_id,
                        !inv ? child_site_id : parent_site_id,
                        partition->repeats->bclv_buffer,
                        inv,
                        partition->attributes);

  free(freqs);
  free(eigenvecs);
  free(inv_eigenvecs);

  return retval;
}

/* computes the table containing the constant parts of the likelihood function
 * partial derivatives on the branch lengths.
 * sumtable: [output] must be allocated for storing (rates x states_padded) values */

PLL_EXPORT int pll_update_sumtable(pll_partition_t * partition,
                                      unsigned int parent_clv_index,
                                      double * parent_clv_vector,
                                      unsigned int child_clv_index,
                                      double * child_clv_vector,
                                      unsigned int * parent_scaler,
                                      unsigned int * child_scaler,
                                      const unsigned int * params_indices,
                                      double *sumtable)
{
  int retval;

  if (pll_repeats_enabled(partition) && 
      (partition->repeats->pernode_ids[parent_clv_index] 
       || partition->repeats->pernode_ids[child_clv_index])) 
  {
    retval = sumtable_repeats(partition,
                                 parent_clv_index,
                                 parent_clv_vector,
                                 child_clv_index,
                                 child_clv_vector,
                                 parent_scaler,
                                 child_scaler,
                                 params_indices,
                                 sumtable);
  }
  else if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
  {
    if ((parent_clv_index < partition->tips) &&
        (child_clv_index < partition->tips))
    {
      /* tip-tip case */
      pll_errno = PLL_ERROR_PARAM_INVALID;
      snprintf(pll_errmsg, 200,
               "pll_update_sumtable() was called for the tip-tip case!");
      retval = PLL_FAILURE;
    }
    else if ((parent_clv_index < partition->tips) ||
             (child_clv_index < partition->tips))
    {
      /* tip-inner */
      retval = sumtable_tipinner(partition,
                                 parent_clv_index,
                                 parent_clv_vector,
                                 child_clv_index,
                                 child_clv_vector,
                                 parent_scaler,
                                 child_scaler,
                                 params_indices,
                                 sumtable);
    }
    else
    {
      /* inner-inner */
      retval = sumtable_innerinner(partition,
                                   parent_clv_index,
                                   parent_clv_vector,
                                   child_clv_index,
                                   child_clv_vector,
                                   parent_scaler,
                                   child_scaler,
                                   params_indices,
                                   sumtable);
    }
  }
  else
  {
    /* inner-inner */
    retval = sumtable_innerinner(partition,
                                 parent_clv_index,
                                 parent_clv_vector,
                                 child_clv_index,
                                 child_clv_vector,
                                 parent_scaler,
                                 child_scaler,
                                 params_indices,
                                 sumtable);
  }

  return retval;
}

PLL_EXPORT void pll_compute_eigenvals_and_prop_invar(pll_partition_t * partition, const unsigned int * params_indices, double *** eigenvals_ptr, double ** prop_invar_ptr)
{
  unsigned int rate_cats = partition->rate_cats;
  double ** eigenvals = (double **) malloc(rate_cats * sizeof(double *));
  double * prop_invar = (double *)  malloc(rate_cats * sizeof(double));
  if (!eigenvals || !prop_invar)
  {
    if (eigenvals) free(eigenvals);
    if (prop_invar) free(prop_invar);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }
  for (unsigned int i=0; i<rate_cats; ++i)
  {
    eigenvals[i]  = partition->eigenvals[params_indices[i]];
    prop_invar[i] = partition->prop_invar[params_indices[i]];
  }
  *eigenvals_ptr = eigenvals;
  *prop_invar_ptr = prop_invar;
}

/* Computes partial derivatives on the branch lengths.
 * branch_length: [input] value where the derivative is computed
 * sumtable: [input] must be computed at the edge where the derivatives will
 *                   be computed
 * f: [output] zeroth derivative
 * d_f:  [output] first derivative
 * dd_f: [output] second derivative
 */
PLL_EXPORT int pll_compute_loglikelihood_derivatives(pll_partition_t * partition,
                                                  int parent_scaler_index,
                                                  unsigned int * parent_scaler,
                                                  int child_scaler_index,
                                                  unsigned int * child_scaler,
                                                  double branch_length,
                                                  const unsigned int * params_indices,
                                                  const double * sumtable,
                                                  double * f,
                                                  double * d_f,
                                                  double * dd_f,
                                                  double * diagptable,
                                                  double * prop_invar)
{
  unsigned int i;
  unsigned int rate_cats = partition->rate_cats;

  double ** freqs     = (double **) malloc(rate_cats * sizeof(double *));
  if (!freqs)
  {
    if (freqs) free(freqs);

    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    return PLL_FAILURE;
  }

  for (i=0; i<rate_cats; ++i)
  {
    freqs[i]      = partition->frequencies[params_indices[i]];
  }

  unsigned int parent_ids = partition->sites;
  unsigned int child_ids = partition->sites;
  if (pll_repeats_enabled(partition))
  {
    parent_ids = parent_scaler_index != PLL_SCALE_BUFFER_NONE 
      ? partition->repeats->perscale_ids[parent_scaler_index]
      : 0;
    parent_ids = parent_ids ? parent_ids : partition->sites;
    child_ids = child_scaler_index != PLL_SCALE_BUFFER_NONE 
      ? partition->repeats->perscale_ids[child_scaler_index]
      : 0;
    child_ids = child_ids ? child_ids : partition->sites;
  }
  int retval = pll_core_loglikelihood_derivatives(partition->states,
                                               partition->sites,
                                               partition->rate_cats,
                                               partition->rate_weights,
                                               parent_scaler,
                                               child_scaler,
                                               parent_ids,
                                               child_ids,
                                               partition->invariant,
                                               partition->pattern_weights,
                                               branch_length,
                                               prop_invar,
                                               freqs,
                                               partition->rates,
                                               sumtable,
                                               f,
                                               d_f,
                                               dd_f,
                                               partition->attributes,
                                               diagptable);

  free (freqs);

  return retval;
}
