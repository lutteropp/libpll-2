/*
 * parse_unetwork_functions.c
 *
 *  Created on: Apr 26, 2019
 *      Author: sarah
 */

#include "pll.h"

/*static pll_unetwork_node_t * alloc_node()
{
  pll_unetwork_node_t * node = (pll_unetwork_node_t *)calloc(1, sizeof(pll_unetwork_node_t));
  node->reticulation_index = -1;
  return node;
}*/

static void dealloc_data(pll_unetwork_node_t * node, void (*cb_destroy)(void *))
{
  if (node->data)
  {
    if (cb_destroy)
      cb_destroy(node->data);
  }
}
/*static void close_roundabout(pll_unetwork_node_t * first)
{
  pll_unetwork_node_t * last = first;
  while(last->next != NULL && last->next != first)
  {
  	if (!last->next->label)
  	  last->next->label = last->label;
  	if (!last->next->reticulation_name)
  	  last->next->reticulation_name = last->reticulation_name;
  	if (last->next->reticulation_index == -1)
  	  last->next->reticulation_index = last->reticulation_index;
  	last = last->next;
  }
  last->next = first;
}*/

static void dealloc_graph_recursive(pll_unetwork_node_t * node,
                                   void (*cb_destroy)(void *),
                                   int level)
{
  if (!node->next)
  {
    /* tip node */
    dealloc_data(node, cb_destroy);
    free(node->label);
    if (node->reticulation_name)
      free(node->reticulation_name);
    free(node);
  }
  else
  {
    /* inner node */
    if (node->label)
      free(node->label);
    if (node->reticulation_name)
      free(node->reticulation_name);

    pll_unetwork_node_t * snode = node;
	do
    {
      if (node != snode || level == 0)
        dealloc_graph_recursive(snode->back, cb_destroy, level+1);
      pll_unetwork_node_t * next = snode->next;
      dealloc_data(snode, cb_destroy);
      free(snode);
      snode = next;
    }
    while(snode && snode != node);
  }
}

PLL_EXPORT void pll_unetwork_graph_destroy(pll_unetwork_node_t * root,
                                        void (*cb_destroy)(void *))
{
  if (!root) return;

  dealloc_graph_recursive(root, cb_destroy, 0);
}

PLL_EXPORT void pll_unetwork_destroy(pll_unetwork_t * network,
                                  void (*cb_destroy)(void *))
{
  unsigned int i;
  pll_unetwork_node_t * node;

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
  if (network->nodes)
  {
    free(network->nodes);
  }
  if (network->reticulation_nodes)
  {
    free(network->reticulation_nodes);
  }
  free(network);
}

static int unetwork_is_rooted(const pll_unetwork_node_t * root)
{
  return (root->next && root->next->next == root) ? 1 : 0;
}

static void fill_nodes_recursive(pll_unetwork_node_t * node,
                                         pll_unetwork_node_t ** array,
		                                 pll_unetwork_node_t ** reticulation_nodes,
                                         unsigned int array_size,
                                         unsigned int * tip_index,
                                         unsigned int * inner_index,
                                         unsigned int level,
	                                     int* visited_reticulations) // TODO: Is this correct?
{
  if (!node || (pll_unetwork_is_reticulation(node) && visited_reticulations[node->reticulation_index])) {
	  return;
  }

  pll_unetwork_node_t * snode = node;
  do
  {
	if (!snode->incoming) // outgoing edge
	{
	  fill_nodes_recursive(snode->back, array, reticulation_nodes, array_size, tip_index, inner_index, level + 1, visited_reticulations);
	}
	snode = snode->next;
  } while (snode && snode != node);

  // now, we are at the node itself.
  unsigned int index;
  if (!node->next)
  {
	index = *tip_index;
	*tip_index += 1;
  }
  else
  { // inner node
	if (pll_unetwork_is_reticulation(node)) {
	  	reticulation_nodes[node->reticulation_index] = node;
	  	visited_reticulations[node->reticulation_index] += 1;
	}
	// TODO: We need to ensure that we only put one representative of each node into the array!!! That is, only the incoming node...
	assert(level == 0 || node->incoming);
	index = *inner_index;
	*inner_index += 1;
  }
  assert(index < array_size);
  array[index] = node;
}

static unsigned int unetwork_count_nodes_recursive(pll_unetwork_node_t * node,
                                                unsigned int * tip_count,
                                                unsigned int * inner_tree_count,
												unsigned int * reticulation_count,
                                                unsigned int level,
												int* visited_reticulations)
{
  if (!node->next) // we have a tip node
  {
    *tip_count += 1;
    return 1;
  }
  else
  {
    unsigned int count = 0;

    pll_unetwork_node_t * snode = level ? node->next : node;
	do
	{
	  if (!snode->incoming)
	  {
		if (!pll_unetwork_is_reticulation(snode->back) || !visited_reticulations[snode->back->reticulation_index])
	      count += unetwork_count_nodes_recursive(snode->back, tip_count, inner_tree_count, reticulation_count, level+1, visited_reticulations);
	  }
	  snode = snode->next;
	}
	while (snode != node);

	if (pll_unetwork_is_reticulation(node)) {
	  *reticulation_count += 1;
	  visited_reticulations[node->reticulation_index] += 1;
	} else {
      *inner_tree_count += 1;
	}

	return count + 1;
  }
}

static unsigned int unetwork_count_nodes(pll_unetwork_node_t * root, unsigned int * tip_count,
                                      unsigned int * inner_tree_count, unsigned int * reticulation_count)
{
  unsigned int count = 0;

  if (tip_count)
    *tip_count = 0;

  if (inner_tree_count)
    *inner_tree_count = 0;

  if (reticulation_count)
	*reticulation_count = 0;

  if (!root->next && !root->back->next)
    return 0;

  if (!root->next)
    root = root->back;

  int* visited_reticulations = (int *)malloc(MAX_RETICULATION_COUNT * sizeof(int));
  int i;
  for (i = 0; i < MAX_RETICULATION_COUNT; ++i)
  {
	visited_reticulations[i] = 0;
  }
  count = unetwork_count_nodes_recursive(root, tip_count, inner_tree_count, reticulation_count, 0, visited_reticulations);
  free(visited_reticulations);

  if (tip_count && inner_tree_count && reticulation_count)
    assert(count == *tip_count + *inner_tree_count + *reticulation_count);

  return count;
}

static void pll_unetwork_unset_indices(pll_unetwork_t * network) {
	unsigned int node_count = network->tip_count + network->inner_tree_count + network->reticulation_count;
	unsigned int i;
	for (i = 0; i < node_count; ++i) {
		pll_unetwork_node_t * snode = network->nodes[i];
		do {
			snode->clv_index = 0;
			snode->node_index = 0;
			snode->scaler_index = 0;
			snode->pmatrix_index = 0;
			snode = snode->next;
		} while (snode && snode != network->nodes[i]);
	}
}

PLL_EXPORT void pll_unetwork_set_indices(pll_unetwork_t * network) {
	unsigned int node_count = network->tip_count + network->inner_tree_count + network->reticulation_count;
	unsigned int i;
	for (i = 0; i < network->tip_count; ++i) {
		network->nodes[i]->clv_index = i;
		network->nodes[i]->node_index = i;
		network->nodes[i]->pmatrix_index = i; // this index is for the edges...
		network->nodes[i]->back->pmatrix_index = i;
		network->nodes[i]->scaler_index = PLL_SCALE_BUFFER_NONE;
	}
	unsigned int clv_idx = network->tip_count;
	unsigned int node_idx = network->tip_count;
	unsigned int pmatrix_idx = network->tip_count;

	for (i = network->tip_count; i < node_count; ++i) {
		pll_unetwork_node_t * snode = network->nodes[i];
		do {
			snode->clv_index = clv_idx;
			snode->node_index = node_idx++;
			snode->scaler_index = snode->clv_index - network->tip_count;
			if (snode->back == network->nodes[0]) {
			} else if (snode->back->pmatrix_index == 0) {
				snode->pmatrix_index = pmatrix_idx;
				snode->back->pmatrix_index = pmatrix_idx++;
				assert(snode->pmatrix_index == snode->back->pmatrix_index);
			} else {
				snode->pmatrix_index = snode->back->pmatrix_index;
			}
			snode = snode->next;
		} while (snode && snode != network->nodes[i]);
		++clv_idx;
	}

	assert(network->vroot->back == network->nodes[0] || network->vroot->pmatrix_index != 0);
	for (i = 0; i < node_count; ++i) {
		assert(network->nodes[i]->pmatrix_index == network->nodes[i]->back->pmatrix_index);
		if (i < network->tip_count) {
			assert(network->nodes[i]->pmatrix_index == network->nodes[i]->clv_index);
			assert(network->nodes[i]->back->pmatrix_index == network->nodes[i]->clv_index);
		}
	}
}

static pll_unetwork_t * unetwork_wrapnetwork(pll_unetwork_node_t * root,
                                    unsigned int tip_count,
                                    unsigned int inner_tree_count,
									unsigned int reticulation_count,
                                    int binary)
{
  unsigned int node_count;

  pll_unetwork_t * network = (pll_unetwork_t *)malloc(sizeof(pll_unetwork_t));
  if (!network)
  {
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    pll_errno = PLL_ERROR_MEM_ALLOC;
    return PLL_FAILURE;
  }

  if (tip_count < 3 && tip_count != 0)
  {
    snprintf(pll_errmsg, 200, "Invalid tip_count value (%u).", tip_count);
    pll_errno = PLL_ERROR_PARAM_INVALID;
    return PLL_FAILURE;
  }

  if (!root->next)
    root = root->back;

  if (binary)
  {
    if (tip_count == 0)
    {
      node_count = unetwork_count_nodes(root, &tip_count, &inner_tree_count, &reticulation_count);
      if (inner_tree_count != tip_count - 2 + reticulation_count)
      {
        snprintf(pll_errmsg, 200, "Input network is not strictly bifurcating.");
        pll_errno = PLL_ERROR_PARAM_INVALID;
        return PLL_FAILURE;
      }
    }
    else
    {
      inner_tree_count = tip_count - 2 + reticulation_count;
      node_count = tip_count + inner_tree_count + reticulation_count;
    }
  }
  else
  {
    if (tip_count == 0 || inner_tree_count == 0)
      node_count = unetwork_count_nodes(root, &tip_count, &inner_tree_count, &reticulation_count);
    else
      node_count = tip_count + inner_tree_count + reticulation_count;
  }

  if (!tip_count)
  {
    snprintf(pll_errmsg, 200, "Input network contains no inner nodes.");
    pll_errno = PLL_ERROR_PARAM_INVALID;
    return PLL_FAILURE;
  }

  network->nodes = (pll_unetwork_node_t **)malloc(node_count*sizeof(pll_unetwork_node_t *));
  if (!network->nodes)
  {
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    pll_errno = PLL_ERROR_MEM_ALLOC;
    return PLL_FAILURE;
  }

  if (reticulation_count > 0)
  {
    network->reticulation_nodes = (pll_unetwork_node_t **)malloc(reticulation_count*sizeof(pll_unetwork_node_t *));
	if (!network->reticulation_nodes)
	{
	  snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
	  pll_errno = PLL_ERROR_MEM_ALLOC;
	  return PLL_FAILURE;
	}
  }
  else
  {
	network->reticulation_nodes = NULL;
  }

  unsigned int tip_index = 0;
  unsigned int inner_index = tip_count;

  int* visited_reticulations = (int *)malloc(MAX_RETICULATION_COUNT * sizeof(int));
  int i;
  for (i = 0; i < MAX_RETICULATION_COUNT; ++i)
  {
  	visited_reticulations[i] = 0;
  }
  fill_nodes_recursive(root, network->nodes, network->reticulation_nodes, node_count, &tip_index, &inner_index, 0, visited_reticulations);
  free(visited_reticulations);

  assert(tip_index == tip_count);
  assert(inner_index == tip_count + inner_tree_count + reticulation_count);

  network->tip_count = tip_count;
  network->inner_tree_count = inner_tree_count;
  network->reticulation_count = reticulation_count;
  network->edge_count = network->tip_count + network->inner_tree_count - 1 + 2 * network->reticulation_count;
  network->tree_edge_count = network->tip_count + network->inner_tree_count - 1; //network->inner_tree_count * 2 + 1; //?
  network->binary = (inner_tree_count == tip_count - (unetwork_is_rooted(root) ? 1 : 2));
  network->vroot = root;

  pll_unetwork_unset_indices(network); // just to be sure
  pll_unetwork_set_indices(network);

  return network;
}

PLL_EXPORT pll_unetwork_t * pll_unetwork_wrapnetwork_multi(pll_unetwork_node_t * root,
                                                  unsigned int tip_count,
                                                  unsigned int inner_tree_count,
												  unsigned int reticulation_count)
{
  return unetwork_wrapnetwork(root, tip_count, inner_tree_count, reticulation_count, 0);
}

/* wraps/encalupsates the unrooted tree graph into a tree structure
   that contains a list of nodes, number of tips and number of inner
   nodes. If 0 is passed as tip_count, then an additional recrursion
   of the tree structure is done to detect the number of tips */
PLL_EXPORT pll_unetwork_t * pll_unetwork_wrapnetwork(pll_unetwork_node_t * root,
                                            unsigned int tip_count)
{
  return unetwork_wrapnetwork(root, tip_count, 0, 0, 1);
}

pll_unetwork_node_t * pll_unetwork_unroot_inplace(pll_unetwork_node_t * root)
{
  /* check for a bifurcation at the root */
  if (unetwork_is_rooted(root))
  {
    pll_unetwork_node_t * left = root->back;
    pll_unetwork_node_t * right =  root->next->back;

    if (root->label)
      free(root->label);
    if (root->reticulation_name)
      free(root->reticulation_name);
    free(root->next);
    free(root);

    double new_length = left->length + right->length;
    left->back = right;
    right->back = left;
    left->length = right->length = new_length;
    left->prob = 1.0;
    left->pmatrix_index = right->pmatrix_index =
        PLL_MIN(left->pmatrix_index, right->pmatrix_index);

    return left->next ? left : right;
  }
  else
  	return root;
}

PLL_EXPORT pll_unetwork_t * pll_unetwork_parse_newick_string(const char * s)
{
  pll_rnetwork_t * rnetwork = pll_rnetwork_parse_newick_string(s);
  assert(rnetwork);
  pll_unetwork_t * unetwork = pll_rnetwork_unroot(rnetwork);
  unetwork->binary = rnetwork->binary;
  pll_rnetwork_destroy(rnetwork, NULL);
  return unetwork;
}

PLL_EXPORT pll_unetwork_t * pll_unetwork_parse_newick(const char * filename)
{
  pll_rnetwork_t * rnetwork = pll_rnetwork_parse_newick(filename);
  pll_unetwork_t * unetwork = pll_rnetwork_unroot(rnetwork);
  pll_rnetwork_destroy(rnetwork, NULL);
  return unetwork;
}
