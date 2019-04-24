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
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "pll.h"

static int unetwork_find(pll_unetwork_node_t * start, pll_unetwork_node_t * target)
{
  /* checks whether the subnetwork rooted at 'start' (in the direction of
     start->next and start->next->next) contains the node 'target' */

  if (!start) return 0;

  if (start == target) return 1;

  if (start->next)
  {
    if (start->next == target) return 1;
    if (unetwork_find(start->next->back, target)) return 1;
  }
  else
    return 0;

  if (start->next->next == target) return 1;
  if (unetwork_find(start->next->next->back, target)) return 1;

  return 0;
}

static void unetwork_link(pll_unetwork_node_t * a,
                       pll_unetwork_node_t * b,
                       double length,
                       unsigned int pmatrix_index)
{
  a->back = b;
  b->back = a;
  a->length = length;
  b->length = length;

  a->pmatrix_index = b->pmatrix_index = pmatrix_index;
}

static void unetwork_swap(pll_unetwork_node_t * t1, pll_unetwork_node_t * t2)
{
  /* swaps the positions of networks t1 and t2. The two networks retain the branch
  lengths from their root to their respective parent nodes, and retain their
  pmatrix indices (i.e. no updating of pmatrices is required) */

  pll_unetwork_node_t * temp = t1->back;

  unetwork_link(t1, t2->back, t2->back->length, t2->back->pmatrix_index);
  unetwork_link(t2, temp, temp->length, temp->pmatrix_index);
}

PLL_EXPORT int pll_unetwork_nni(pll_unetwork_node_t * p,
                             int type,
                             pll_unetwork_rb_t * rb)
{
  pll_unetwork_node_t * subnetwork1;
  pll_unetwork_node_t * subnetwork2;

  if ((type != PLL_UNETWORK_MOVE_NNI_LEFT) && (type != PLL_UNETWORK_MOVE_NNI_RIGHT))
  {
    snprintf(pll_errmsg, 200, "Invalid NNI move type");
    pll_errno = PLL_ERROR_NNI_INVALIDMOVE;
    return PLL_FAILURE;
  }

  /* check if selected node p is edge  */
  if (!(p->next) || !(p->back->next))
  {
    snprintf(pll_errmsg, 200, "Specified terminal branch");
    pll_errno = PLL_ERROR_NNI_TERMINALBRANCH;
    return PLL_FAILURE;
  }

  /* check if rollback buffer is provided, and fill it up */
  if (rb)
  {
    rb->move_type = PLL_UNETWORK_MOVE_NNI;
    rb->nni.p = p;
    rb->nni.nni_type = type;
  }

  subnetwork1 = p->next;
  subnetwork2 = (type == PLL_UNETWORK_MOVE_NNI_LEFT) ?
               p->back->next : p->back->next->next;

  unetwork_swap(subnetwork1,subnetwork2);

  return PLL_SUCCESS;
}

static int unetwork_nni_rollback(pll_unetwork_rb_t * rb)
{
  /* restore the network topology from a previous SPR */
  return pll_unetwork_nni(rb->nni.p,
                       rb->nni.nni_type,
                       NULL);
}

PLL_EXPORT int pll_unetwork_spr(pll_unetwork_node_t * p,
                             pll_unetwork_node_t * r,
                             pll_unetwork_rb_t * rb,
                             double * branch_lengths,
                             unsigned int * matrix_indices)
{
  /* given nodes p and r, perform an SPR move in the following way,
     i.e. prune subnetwork C and make it adjacent to subnetwork D:

      A           B          C             D           A          B
     ____        ____       ____          ____        ____       ____
     \  /        \  /       \  /          \  /        \  /       \  /
      \/          \/         \/            \/          \/         \/
       *          *          * p'           *          *          *
        \         |     q   /                \         |         /
         *'*_____.*._____*'* p     --->       *'*_____.*._____*'*
         '*'     *.*     '*'                  '*'     *.*     '*'
         / r       u    q' \                  /                 \
     r' *                   * v              *                   *
       /\                   /\              /\                   /\
      /__\                 /__\            /__\                 /__\

       D                    E               C                    E

     node p must be part of an inner node (i.e. node with ->next set). The
     procedure prunes the subnetwork rooted at the opposite end-point of p
     (subnetwork C in our case) and regrafts it on the edge r'<->r. It is done
     in the following way:

     (a) prune the subnetwork rooted at the opposite end-point of p (p' on figure)
         by breaking the edges q<->u and q'<->v

     (b) connect node u with node v

     (c) break edge r<->r' by connecting node r with node q, and node r' with
         node q'

     Node r must not be part of the subnetwork to be pruned (C in this case). Note
     that for speed reasons, the function *does not* check this property to save
     a network traversal. A safer (albeit slower) function that checks this
     property is pll_unetwork_spr_safe
  */

  int k = 0;

  if ((!branch_lengths && matrix_indices) ||
      (branch_lengths && !matrix_indices))
  {
    pll_errno = PLL_ERROR_PARAM_INVALID;
    snprintf(pll_errmsg, 200, "Parameters 4,5 must be both NULL or both set");
    return PLL_FAILURE;
  }

  /* if p is a tip node then prompt an error */
  if (!p->next)
  {
    pll_errno = PLL_ERROR_SPR_TERMINALBRANCH;
    snprintf(pll_errmsg, 200, "Prune edge must be defined by an inner node");
    return PLL_FAILURE;
  }

  /* check whether the move will result in the same network */
  if (r == p || r == p->back ||
      r == p->next || r == p->next->back ||
      r == p->next->next || r == p->next->next->back)
  {
    pll_errno = PLL_ERROR_SPR_NOCHANGE;
    snprintf(pll_errmsg, 200, "Proposed move yields the same network");
    return PLL_FAILURE;
  }

  /* check if rollback buffer is provided, and fill it up */
  if (rb)
  {
    rb->move_type = PLL_UNETWORK_MOVE_SPR;
    rb->spr.p = p;
    rb->spr.r = r;
    rb->spr.rb = r->back;
    rb->spr.r_len = r->length;
    rb->spr.pnb = p->next->back;
    rb->spr.pnb_len = p->next->length;
    rb->spr.pnnb = p->next->next->back;
    rb->spr.pnnb_len = p->next->next->length;
  }

  /* (b) connect u and v */
  pll_unetwork_node_t * u = p->next->back;
  pll_unetwork_node_t * v = p->next->next->back;
  unetwork_link(u,
             v,
             u->length + v->length,
             u->pmatrix_index);
  /* if requested, store the new branch length for the corresponding 
     pmatrix index */
  if (branch_lengths)
  {
    branch_lengths[k] = u->length;
    matrix_indices[k] = u->pmatrix_index;
  }

  /* (a) prune subnetwork C */
  p->next->back = p->next->next->back = NULL;

  /* (c) regraft C at r<->r' */
  double length = r->length / 2;

  /* r' <-> q' */
  unetwork_link(r->back,
             p->next->next,
             length,
             p->next->next->pmatrix_index);
  /* if requested, store the new branch length for the corresponding 
     pmatrix index */
  if (branch_lengths)
  {
    ++k;
    branch_lengths[k] = length;
    matrix_indices[k]   = p->next->next->pmatrix_index;
  }

  /* r<->q */
  unetwork_link(r,
             p->next,
             length,
             r->pmatrix_index);
  /* if requested, store the new branch length for the corresponding 
     pmatrix index */
  if (branch_lengths)
  {
    ++k;
    branch_lengths[k] = length;
    matrix_indices[k]   = r->pmatrix_index;
  }

  return PLL_SUCCESS;
}

static int unetwork_spr_rollback(pll_unetwork_rb_t * rb,
                              double * branch_lengths,
                              unsigned int * matrix_indices)
{
  if ((!branch_lengths && matrix_indices) ||
      (branch_lengths && !matrix_indices))
  {
    pll_errno = PLL_ERROR_PARAM_INVALID;
    snprintf(pll_errmsg, 200, "Parameters 4,5 must be both NULL or both set");
    return PLL_FAILURE;
  }

  int k = 0;

  /* restore the network topology from a previous SPR */
  unetwork_link(rb->spr.pnb,
             rb->spr.p->next,
             rb->spr.pnb_len,
             rb->spr.pnb->pmatrix_index);
  if (branch_lengths)
  {
    branch_lengths[k] = rb->spr.pnb_len;
    matrix_indices[k] = rb->spr.pnb->pmatrix_index;
  }

  unetwork_link(rb->spr.pnnb,
             rb->spr.p->next->next,
             rb->spr.pnnb_len,
             rb->spr.p->next->next->pmatrix_index);
  if (branch_lengths)
  {
    branch_lengths[++k] = rb->spr.pnnb_len;
    matrix_indices[k]   = rb->spr.p->next->next->pmatrix_index;
  }

  unetwork_link(rb->spr.r,
             rb->spr.rb,
             rb->spr.r_len,
             rb->spr.r->pmatrix_index);
  if (branch_lengths)
  {
    branch_lengths[++k] = rb->spr.r_len;
    matrix_indices[k]   = rb->spr.r->pmatrix_index;
  }

  return PLL_SUCCESS;
}

/* this is a safer (but slower) function for performing an spr move, than
   pll_unetwork_spr(). See the last paragraph in the comments section of the
   pll_unetwork_spr() function for more details */
PLL_EXPORT int pll_unetwork_spr_safe(pll_unetwork_node_t * p,
                                  pll_unetwork_node_t * r,
                                  pll_unetwork_rb_t * rb,
                                  double * branch_lengths,
                                  unsigned int * matrix_indices)
{
  /* check all possible scenarios of failure */
  if (!p)
  {
    pll_errno = PLL_ERROR_PARAM_INVALID;
    snprintf(pll_errmsg, 200, "Node p is set to NULL");
    return PLL_FAILURE;
  }

  if (!r)
  {
    pll_errno = PLL_ERROR_PARAM_INVALID;
    snprintf(pll_errmsg, 200, "Node r is set to NULL");
    return PLL_FAILURE;
  }

  if (!p->next)
  {
    pll_errno = PLL_ERROR_SPR_TERMINALBRANCH;
    snprintf(pll_errmsg, 200, "Prune edge must be defined by an inner node");
    return PLL_FAILURE;
  }

  /* check whether the move results in the same network */
  if (r == p || r == p->back ||
      r == p->next || r == p->next->back ||
      r == p->next->next || r == p->next->next->back)
  {
    pll_errno = PLL_ERROR_SPR_NOCHANGE;
    snprintf(pll_errmsg, 200, "Proposed move yields the same network");
    return PLL_FAILURE;
  }

  /* node r must not be in the same subnetwork as the one that is to be pruned */
  if (unetwork_find(p->back, r))
  {
    pll_errno = PLL_ERROR_PARAM_INVALID;
    snprintf(pll_errmsg, 200, "Node r is part of the subnetwork to be pruned");
    return PLL_FAILURE;
  }

  return pll_unetwork_spr(p,r,rb,branch_lengths,matrix_indices);
}

PLL_EXPORT int pll_unetwork_rollback(pll_unetwork_rb_t * rollback,
                                  double * branch_lengths,
                                  unsigned int * matrix_indices)
{
  if (!rollback)
  {
    pll_errno = PLL_ERROR_PARAM_INVALID;
    snprintf(pll_errmsg, 200, "Provide a rollback");
    return PLL_FAILURE;
  }

  if (rollback->move_type == PLL_UNETWORK_MOVE_SPR)
    return unetwork_spr_rollback(rollback, branch_lengths, matrix_indices);
  else if (rollback->move_type == PLL_UNETWORK_MOVE_NNI)
    return unetwork_nni_rollback(rollback);

  pll_errno = PLL_ERROR_PARAM_INVALID;
  snprintf(pll_errmsg, 200, "Invalid move type");
  return PLL_FAILURE;
}
