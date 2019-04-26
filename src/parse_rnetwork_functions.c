/*
 * parse_rnetwork_functions.c
 *
 *  Created on: Apr 26, 2019
 *      Author: sarah
 */

#include "pll.h"

/*static void recursive_assign_indices(pll_rnetwork_node_t * node,
                                     unsigned int * tip_clv_index,
                                     unsigned int * inner_clv_index,
                                     int * inner_scaler_index,
                                     unsigned int * inner_node_index)
{
  if (!node->left)
  {
    node->node_index = *tip_clv_index;
    node->clv_index = *tip_clv_index;
    node->pmatrix_index = *tip_clv_index;
    node->scaler_index = PLL_SCALE_BUFFER_NONE;
    *tip_clv_index = *tip_clv_index + 1;
    return;
  }

  recursive_assign_indices(node->left,
                           tip_clv_index,
                           inner_clv_index,
                           inner_scaler_index,
                           inner_node_index);

  recursive_assign_indices(node->right,
                           tip_clv_index,
                           inner_clv_index,
                           inner_scaler_index,
                           inner_node_index);

  node->node_index = *inner_node_index;
  node->clv_index = *inner_clv_index;
  node->scaler_index = *inner_scaler_index;
  node->pmatrix_index = *inner_clv_index;

  *inner_clv_index = *inner_clv_index + 1;
  *inner_scaler_index = *inner_scaler_index + 1;
  *inner_node_index = *inner_node_index + 1;
}*/

/*
PLL_EXPORT void pll_rnetwork_reset_template_indices(pll_rnetwork_node_t * root,
                                                 unsigned int tip_count)
{
  unsigned int tip_clv_index = 0;
  unsigned int inner_clv_index = tip_count;
  unsigned int inner_node_index = tip_count;
  int inner_scaler_index = 0;

  recursive_assign_indices(root->left,
                           &tip_clv_index,
                           &inner_clv_index,
                           &inner_scaler_index,
                           &inner_node_index);

  recursive_assign_indices(root->right,
                           &tip_clv_index,
                           &inner_clv_index,
                           &inner_scaler_index,
                           &inner_node_index);

  root->node_index = inner_node_index;
  root->clv_index = inner_clv_index;
  root->scaler_index = inner_scaler_index;

  // root gets any number for pmatrix since it will never be used
  root->pmatrix_index = 0;
}
*/

static void fill_nodes_recursive(pll_rnetwork_node_t * node,
                                 pll_rnetwork_node_t ** array,
                                 pll_rnetwork_node_t ** reticulations,
                                 unsigned int * tip_index,
                                 unsigned int * inner_index,
                                 unsigned int * scaler_index)
{
  if (!node)
  {
    return;
  }

  if (!node->is_reticulation && !node->left && !node->right)
  { // we are at a tip node
    node->idx = *tip_index;
    node->clv_index = *tip_index;
    node->pmatrix_index = *tip_index;
    node->scaler_index = PLL_SCALE_BUFFER_NONE;
    array[*tip_index] = node;
    *tip_index = *tip_index + 1;
    return;
  }

  if (node->idx > 0) // the node has already been visited
  {
    return;
  }

  if (!node->is_reticulation)
  {
    fill_nodes_recursive(node->left,  array, reticulations, tip_index, inner_index, scaler_index);
    fill_nodes_recursive(node->right, array, reticulations, tip_index, inner_index, scaler_index);
  }
  else
  {
    reticulations[node->reticulation_index] = node;
    fill_nodes_recursive(node->child, array, reticulations, tip_index, inner_index, scaler_index);
  }

  array[*inner_index] = node;
  node->idx = *inner_index;
  node->clv_index = *inner_index;
  node->pmatrix_index = *inner_index;
  node->scaler_index = *scaler_index;
  *inner_index = *inner_index + 1;
  *scaler_index = *scaler_index + 1;
}

static unsigned int rnetwork_count_nodes_recursive(pll_rnetwork_node_t * node,
                                                unsigned int * tip_count,
                                                unsigned int * inner_tree_count,
												unsigned int * reticulation_count)
{
  if (!node)
  {
	return 0;
  }

  if (!node->is_reticulation)
  {
	 if (!node->left && !node->right)
	 {
	   *tip_count += 1;
	   return 1;
	 }
	 else
	 {
	   *inner_tree_count += 1;
	   return 1 + rnetwork_count_nodes_recursive(node->left, tip_count, inner_tree_count, reticulation_count) + rnetwork_count_nodes_recursive(node->right, tip_count, inner_tree_count, reticulation_count);
	 }
  }
  else
  {
	*reticulation_count += 1;
	return 1 + rnetwork_count_nodes_recursive(node->child, tip_count, inner_tree_count, reticulation_count);
  }
}

static unsigned int rnetwork_count_nodes(pll_rnetwork_node_t * root, unsigned int * tip_count,
                                      unsigned int * inner_tree_count, unsigned int * reticulation_count)
{
  unsigned int count = 0;

  if (tip_count)
    *tip_count = 0;

  if (inner_tree_count)
    *inner_tree_count = 0;

  if (reticulation_count)
	*reticulation_count = 0;

  if (!root)
	return 0;

  count = rnetwork_count_nodes_recursive(root, tip_count, inner_tree_count, reticulation_count);

  if (tip_count && inner_tree_count && reticulation_count)
    assert(count == *tip_count + *inner_tree_count + *reticulation_count);

  return count;
}

PLL_EXPORT pll_rnetwork_t * pll_rnetwork_wrapnetwork(pll_rnetwork_node_t * root)
{
  pll_rnetwork_t * network = (pll_rnetwork_t *)malloc(sizeof(pll_rnetwork_t));
  if (!network)
  {
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    pll_errno = PLL_ERROR_MEM_ALLOC;
    return PLL_FAILURE;
  }

  unsigned int tip_count = 0;
  unsigned int inner_tree_count = 0;
  unsigned int reticulation_count = 0;
  unsigned int node_count = rnetwork_count_nodes(root, &tip_count, &inner_tree_count, &reticulation_count);
  network->tip_count = tip_count;
  network->inner_tree_count = inner_tree_count;
  network->reticulation_count = reticulation_count;

  network->nodes = (pll_rnetwork_node_t **)malloc(node_count * sizeof(pll_rnetwork_node_t *));
  network->reticulation_nodes = (pll_rnetwork_node_t **)malloc(reticulation_count * sizeof(pll_rnetwork_node_t *));
  if (!network->nodes || !network->reticulation_nodes)
  {
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    pll_errno = PLL_ERROR_MEM_ALLOC;
    return PLL_FAILURE;
  }

  unsigned int tip_index = 0;
  unsigned int inner_index = tip_count;
  unsigned int scaler_index = 0;

  fill_nodes_recursive(root->left, network->nodes, network->reticulation_nodes, &tip_index, &inner_index, &scaler_index);
  fill_nodes_recursive(root->right, network->nodes, network->reticulation_nodes, &tip_index, &inner_index, &scaler_index);
  root->idx = inner_index;
  root->clv_index = inner_index;
  root->pmatrix_index = inner_index;
  root->scaler_index = scaler_index;
  network->nodes[inner_index] = root;
  network->root = root;
  network->edge_count = reticulation_count + inner_tree_count * 2;
  network->tree_edge_count = inner_tree_count * 2;
  network->binary = 1;

  return network;
}

static void dealloc_data(pll_rnetwork_node_t * node, void (*cb_destroy)(void *))
{
  if (node->data)
  {
    if (cb_destroy)
      cb_destroy(node->data);
  }
}

PLL_EXPORT void pll_rnetwork_graph_destroy(pll_rnetwork_node_t * root,
                                        void (*cb_destroy)(void *))
{
  if (!root) return;

  pll_rnetwork_graph_destroy(root->left, cb_destroy);
  pll_rnetwork_graph_destroy(root->right, cb_destroy);

  dealloc_data(root, cb_destroy);
  free(root->label);
  free(root);
}

PLL_EXPORT void pll_rnetwork_destroy(pll_rnetwork_t * network,
                                  void (*cb_destroy)(void *))
{
  unsigned int i;
  pll_rnetwork_node_t * node;

  /* deallocate all nodes */
  for (i = 0; i < network->tip_count + network->inner_tree_count + network->reticulation_count; ++i)
  {
    node = network->nodes[i];
    dealloc_data(node, cb_destroy);

    if (node->label)
      free(node->label);
    if (node->reticulation_name)
      free(node->reticulation_name);

    free(node);
  }

  /* deallocate network structure */
  free(network->nodes);
  free(network->reticulation_nodes);
  free(network);
}
