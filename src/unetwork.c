/*
 * unetwork.c
 *
 *  Created on: Apr 23, 2019
 *      Author: sarah
 */

#include "pll.h"
#include <stdarg.h>

__attribute__((format(printf, 1, 2)))
void fatal(const char * format, ...)
{
  va_list argptr;
  va_start(argptr, format);
  vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}

void * xmalloc(size_t size)
{
  void * t;
  t = malloc(size);
  if (!t)
    fatal("Unable to allocate enough memory.");

  return t;
}

char * xstrdup(const char * s)
{
  size_t len = strlen(s);
  char * p = (char *)xmalloc(len+1);
  return strcpy(p,s);
}

unsigned int count_outgoing(const pll_unetwork_node_t * node) {
	unsigned int cnt = 0;
	pll_unetwork_node_t * snode = node->next;
	do {
		if (!snode->incoming) {
			cnt++;
		}
		snode = snode->next;
	} while (snode && snode != node);
	return cnt;
}

unsigned int count_active_outgoing(const pll_unetwork_node_t * node) {
	unsigned int cnt = 0;
	pll_unetwork_node_t * snode = node->next;
	do {
		if (!snode->incoming && snode->active) {
			cnt++;
		}
		snode = snode->next;
	} while (snode && snode != node);
	return cnt;
}

unsigned int count_incoming(const pll_unetwork_node_t * node) {
	unsigned int cnt = 0;
	pll_unetwork_node_t * snode = node->next;
	do {
		if (snode->incoming) {
			cnt++;
		}
		snode = snode->next;
	} while (snode && snode != node);
	return cnt;
}

unsigned int count_active_incoming(const pll_unetwork_node_t * node) {
	unsigned int cnt = 0;
	pll_unetwork_node_t * snode = node->next;
	do {
		if (snode->incoming && snode->active) {
			cnt++;
		}
		snode = snode->next;
	} while (snode && snode != node);
	return cnt;
}

PLL_EXPORT int node_is_inner_tree(const pll_unetwork_node_t * node) {
	unsigned int cnt_out = count_outgoing(node);
	return (cnt_out > 1);
}

PLL_EXPORT int node_is_reticulation(const pll_unetwork_node_t * node) {
	unsigned int cnt_in = count_incoming(node);
	return (cnt_in > 1);
}

PLL_EXPORT int node_is_leaf(const pll_unetwork_node_t * node) {
	unsigned int cnt_out = count_outgoing(node);
	return (cnt_out == 0);
}

PLL_EXPORT int node_is_root(const pll_unetwork_node_t * node) {
	unsigned int cnt_in = count_incoming(node);
	return (cnt_in == 0);
}

PLL_EXPORT char * pll_unetwork_export_newick(const pll_unetwork_node_t * root,
                                   char * (*cb_serialize)(const pll_unetwork_node_t *)) {
	return PLL_FAILURE;
}

PLL_EXPORT char * pll_unetwork_export_newick_rooted(const pll_unetwork_node_t * root,
                                                 double root_brlen) {
	return PLL_FAILURE;
}

static void unetwork_tree_traverse_recursive(pll_unetwork_node_t * node,
                                     int traversal,
                                     int (*cbtrav)(pll_unetwork_node_t *),
                                     unsigned int * index,
									 pll_unetwork_node_t ** outbuffer)
{
  if (!cbtrav(node))
    return;

  if (traversal == PLL_TREE_TRAVERSE_PREORDER)
  {
	if (count_active_outgoing(node) > 1) {
      outbuffer[*index] = node;
      *index = *index + 1;
	}
  }

  if (node->next)
  {
	pll_unetwork_node_t * snode = node->next;
    do
    {
      if (snode->active) {
        unetwork_tree_traverse_recursive(snode->back, traversal, cbtrav, index, outbuffer);
      }
      snode = snode->next;
    }
    while (snode && snode != node);
  }

  if (traversal == PLL_TREE_TRAVERSE_POSTORDER)
  {
	if (count_active_outgoing(node) > 1) {
      outbuffer[*index] = node;
      *index = *index + 1;
	}
  }
}

PLL_EXPORT int pll_unetwork_set_reticulation_parents(pll_unetwork_t * network, uint64_t tree_number) {
  if (!network->binary) {
    return PLL_FAILURE;
  }
  if (tree_number >= ((unsigned int) 2 << network->reticulation_count))
  	return PLL_FAILURE;

  unsigned int i;
  for (i = 0; i < network->reticulation_count; ++i) {
    int take_first_parent = (tree_number >> network->reticulation_nodes[i]->reticulation_index) & 1;
    pll_unetwork_node_t * snode = network->reticulation_nodes[i]->next;
    pll_unetwork_node_t * first_incoming = NULL;
    pll_unetwork_node_t * second_incoming = NULL;
    while (snode && snode != network->reticulation_nodes[i]) {
      if (snode->incoming) {
    	if (!first_incoming) {
    	  first_incoming = snode;
    	} else {
    	  second_incoming = snode;
    	}
      }
      snode = snode->next;
    }
    assert(first_incoming && second_incoming);
    if (take_first_parent) {
      first_incoming->active = 1;
      first_incoming->back->active = 1;
      second_incoming->active = 0;
      second_incoming->back->active = 0;
    } else {
      first_incoming->active = 0;
      first_incoming->back->active = 0;
      second_incoming->active = 1;
      second_incoming->back->active = 1;
    }
  }
  return PLL_SUCCESS;
}
PLL_EXPORT int pll_unetwork_forget_reticulation_parents(pll_unetwork_t * network) {
  unsigned int i;
  for (i = 0; i < network->reticulation_count; ++i) {
	pll_unetwork_node_t * snode = network->reticulation_nodes[i]->next;
	while (snode && snode != network->reticulation_nodes[i]) {
		snode->active = 1;
		snode->back->active = 1;
		snode = snode->next;
	}
  }
  return PLL_SUCCESS;
}

PLL_EXPORT int pll_unetwork_tree_traverse(pll_unetwork_t * network,
                                  int traversal,
                                  int (*cbtrav)(pll_unetwork_node_t *),
                                  pll_unetwork_node_t ** outbuffer,
                                  unsigned int * trav_size, uint64_t tree_number) {
  if (tree_number >= ((unsigned int) 2 << network->reticulation_count))
			return PLL_FAILURE;
  pll_unetwork_set_reticulation_parents(network, tree_number);

  *trav_size = 0;
  if (!network->vroot->next) return PLL_FAILURE;

  if (traversal == PLL_TREE_TRAVERSE_POSTORDER ||
	  traversal == PLL_TREE_TRAVERSE_PREORDER)
  {

	/* we will traverse an unrooted network in the following way

				2
			  /
		1  --*
			  \
				3

	   at each node the callback function is called to decide whether we
	   are going to traversing the subtree rooted at the specific node */

	if (network->vroot->back->active) {
	  unetwork_tree_traverse_recursive(network->vroot->back, traversal, cbtrav, trav_size, outbuffer);
	}
	if (network->vroot->active) {
	  unetwork_tree_traverse_recursive(network->vroot, traversal, cbtrav, trav_size, outbuffer);
	}
  }
  else
  {
	snprintf(pll_errmsg, 200, "Invalid traversal value.");
	pll_errno = PLL_ERROR_PARAM_INVALID;
	return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}

PLL_EXPORT void pll_unetwork_create_operations(pll_unetwork_node_t * const* trav_buffer,
                                            unsigned int trav_buffer_size,
                                            double * branches,
                                            unsigned int * pmatrix_indices,
                                            pll_operation_t * ops,
                                            unsigned int * matrix_count,
                                            unsigned int * ops_count) {
	return;
}

/* a callback function for checking network tree integrity */
static int cb_check_tree_integrity_mult(const pll_unetwork_t * network,
                                   const pll_unetwork_node_t * node)
{
  unsigned int clv_index = node->clv_index;
  int scaler_index = node->scaler_index;
  unsigned int pmatrix_index = node->pmatrix_index;
  char * label = node->label;
  double length = node->length;
  unsigned int subnodes = 1;

  if (count_active_incoming(node) > 1) {
	  snprintf(pll_errmsg, 200, "Encountered a node with more than one active incoming edge");
	  return PLL_FAILURE;
  }

  /* edge attributes */
  if (node->back->length != length)
  {
    snprintf(pll_errmsg, 200, "Inconsistent branch lengths: %lf != %lf",
             length, node->back->length);
    return PLL_FAILURE;
  }

  if (node->back->pmatrix_index != pmatrix_index)
  {
    snprintf(pll_errmsg, 200, "Inconsistent pmatrix indices: %u != %u",
             pmatrix_index, node->back->pmatrix_index);
    return PLL_FAILURE;
  }

  if (node->next)
  {
    /* node attributes */
    pll_unetwork_node_t * snode = node->next;
    do
    {
      subnodes++;

      if (network->binary && subnodes > 3 && !node_is_reticulation(snode))
      {
        snprintf(pll_errmsg, 200, "Multifurcation found in a binary network "
            "at node with clv_index = %u", snode->clv_index);
        return PLL_FAILURE;
      }

      if (subnodes > network->tip_count)
      {
        snprintf(pll_errmsg, 200, "Multifurcation exceeding the network size found "
            "at node with clv_index = %u", snode->clv_index);
        return PLL_FAILURE;
      }

      if (snode->clv_index != clv_index)
      {
        snprintf(pll_errmsg, 200, "Inconsistent CLV indices: %u != %u",
                 clv_index, snode->clv_index);
        return PLL_FAILURE;
      }
      if (snode->scaler_index != scaler_index)
      {
        snprintf(pll_errmsg, 200, "Inconsistent scaler indices: %u != %u",
                 scaler_index, snode->scaler_index);
        return PLL_FAILURE;
      }
      if (snode->label != label)
      {
        snprintf(pll_errmsg, 200, "Inconsistent node labels: '%s' != '%s'",
                 label, snode->label);
        return PLL_FAILURE;
      }
      if (!snode->next)
      {
        snprintf(pll_errmsg, 200, "Open roundabout (node->next is NULL) "
            "at node with clv_index = %u", snode->clv_index);
        return PLL_FAILURE;
      }
      snode = snode->next;
    }
    while (snode != node);
  }

  return 1;
}

/* a callback function for checking network integrity */
static int cb_check_integrity_mult(const pll_unetwork_t * network,
                                   const pll_unetwork_node_t * node)
{
  unsigned int clv_index = node->clv_index;
  int scaler_index = node->scaler_index;
  unsigned int pmatrix_index = node->pmatrix_index;
  char * label = node->label;
  double length = node->length;
  unsigned int subnodes = 1;

  /* edge attributes */
  if (node->back->length != length)
  {
    snprintf(pll_errmsg, 200, "Inconsistent branch lengths: %lf != %lf",
             length, node->back->length);
    return PLL_FAILURE;
  }

  if (node->back->pmatrix_index != pmatrix_index)
  {
    snprintf(pll_errmsg, 200, "Inconsistent pmatrix indices: %u != %u",
             pmatrix_index, node->back->pmatrix_index);
    return PLL_FAILURE;
  }

  if (node->next)
  {
    /* node attributes */
    pll_unetwork_node_t * snode = node->next;
    do
    {
      subnodes++;

      if (network->binary && subnodes > 3 && !node_is_reticulation(snode))
      {
        snprintf(pll_errmsg, 200, "Multifurcation found in a binary network "
            "at node with clv_index = %u", snode->clv_index);
        return PLL_FAILURE;
      }

      if (subnodes > network->tip_count)
      {
        snprintf(pll_errmsg, 200, "Multifurcation exceeding the network size found "
            "at node with clv_index = %u", snode->clv_index);
        return PLL_FAILURE;
      }

      if (snode->clv_index != clv_index)
      {
        snprintf(pll_errmsg, 200, "Inconsistent CLV indices: %u != %u",
                 clv_index, snode->clv_index);
        return PLL_FAILURE;
      }
      if (snode->scaler_index != scaler_index)
      {
        snprintf(pll_errmsg, 200, "Inconsistent scaler indices: %u != %u",
                 scaler_index, snode->scaler_index);
        return PLL_FAILURE;
      }
      if (snode->label != label)
      {
        snprintf(pll_errmsg, 200, "Inconsistent node labels: '%s' != '%s'",
                 label, snode->label);
        return PLL_FAILURE;
      }
      if (!snode->next)
      {
        snprintf(pll_errmsg, 200, "Open roundabout (node->next is NULL) "
            "at node with clv_index = %u", snode->clv_index);
        return PLL_FAILURE;
      }
      snode = snode->next;
    }
    while (snode != node);
  }

  return 1;
}

PLL_EXPORT int pll_unetwork_every(pll_unetwork_t * network,
                               int (*cb)(const pll_unetwork_t *,
                                         const pll_unetwork_node_t *))
{
  unsigned int i;
  int rc = 1;

  for (i = 0; i < network->tip_count + network->inner_tree_count + network->reticulation_count; ++i)
    rc &= cb(network, network->nodes[i]);

  return (rc ? PLL_SUCCESS : PLL_FAILURE);
}

PLL_EXPORT int pll_unetwork_every_const(const pll_unetwork_t * network,
                               int (*cb)(const pll_unetwork_t *,
                                         const pll_unetwork_node_t *))
{
  unsigned int i;
  int rc = 1;

  for (i = 0; i < network->tip_count + network->inner_tree_count + network->reticulation_count; ++i)
    rc &= cb(network, network->nodes[i]);

  return (rc ? PLL_SUCCESS : PLL_FAILURE);
}

PLL_EXPORT int pll_unetwork_check_integrity(const pll_unetwork_t * network) {
  return pll_unetwork_every_const(network, cb_check_integrity_mult);
}

PLL_EXPORT int pll_unetwork_check_tree_integrity(const pll_unetwork_t * network) {
  return pll_unetwork_every_const(network, cb_check_tree_integrity_mult);
}

/* TODO: Memory allocation checks were not implemented in this function!!! */
static pll_unetwork_node_t * clone_node(const pll_unetwork_node_t * node)
{
	pll_unetwork_node_t * new_node = (pll_unetwork_node_t *)malloc(sizeof(pll_unetwork_node_t));
  memcpy(new_node, node, sizeof(pll_unetwork_node_t));

  if (node->label)
  {
    new_node->label = (char *)malloc(strlen(node->label)+1);
    strcpy(new_node->label,node->label);
  }
  if (node->reticulation_name)
  {
	new_node->reticulation_name = (char *)malloc(strlen(node->reticulation_name)+1);
	strcpy(new_node->reticulation_name, node->reticulation_name);
  }

  if (node->next)
  {
	pll_unetwork_node_t * snode = node->next;
	pll_unetwork_node_t * new_snode = new_node;
    do
    {
      new_snode->next = (pll_unetwork_node_t *)malloc(sizeof(pll_unetwork_node_t));
      memcpy(new_snode->next, snode, sizeof(pll_unetwork_node_t));
      new_snode->next->label = new_node->label;
      new_snode->next->reticulation_name = new_node->reticulation_name;
      snode = snode->next;
      new_snode = new_snode->next;
    }
    while (snode != node);

    new_snode->next = new_node;
  }

  return new_node;
}

static void unetwork_recurse_clone(pll_unetwork_node_t * new_root, const pll_unetwork_node_t * root)
{
  const pll_unetwork_node_t * node = root->back;
  if (node)
  {
    new_root->back = clone_node(node);
    new_root->back->back = new_root;

    if (node->next)
    {
      pll_unetwork_node_t * snode = node->next;
      pll_unetwork_node_t * new_snode = new_root->back->next;
      do
      {
        unetwork_recurse_clone(new_snode, snode);
        snode = snode->next;
        new_snode = new_snode->next;
      }
      while (snode && snode != node);
    }
  }
}

PLL_EXPORT pll_unetwork_node_t * pll_unetwork_graph_clone(const pll_unetwork_node_t * root) {
  pll_unetwork_node_t * new_root = clone_node(root);

  const pll_unetwork_node_t * snode = root;
  pll_unetwork_node_t * new_snode = new_root;
  do
  {
    unetwork_recurse_clone(new_snode, snode);
	snode = snode->next;
	new_snode = new_snode->next;
  }
  while (snode && snode != root);

  return new_root;
}

PLL_EXPORT pll_unetwork_t * pll_unetwork_clone(const pll_unetwork_t * network) {
   /* choose the last inner node as the starting point of the clone. It does not
	 really matter which node to choose, but since the newick parser places the
	 root node at the end of the list, we use the same notation here */
  pll_unetwork_node_t * root = pll_unetwork_graph_clone(network->vroot);

  if (network->binary)
	return pll_unetwork_wrapnetwork(root, network->tip_count);
  else
	return pll_unetwork_wrapnetwork_multi(root, network->tip_count, network->inner_tree_count, network->reticulation_count);
}

static pll_unetwork_node_t * rnetwork_unroot(pll_rnetwork_node_t * root, pll_unetwork_node_t * back, pll_unetwork_node_t** reticulation_nodes)
{
  pll_unetwork_node_t * uroot;
  if (!root->is_reticulation) {
	  uroot = (void *)calloc(1,sizeof(pll_unetwork_node_t));
	  if (!uroot)
	  {
	    pll_errno = PLL_ERROR_MEM_ALLOC;
	    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
	    return NULL;
	  }
	  uroot->back = back;
	  uroot->label = (root->label) ? xstrdup(root->label) : NULL;
	  uroot->reticulation_name = (root->reticulation_name) ? xstrdup(root->reticulation_name) : NULL;
	  uroot->reticulation_index = -1;
	  uroot->clv_index = root->clv_index;
	  uroot->node_index = root->idx;
	  uroot->scaler_index = root->scaler_index;
	  uroot->pmatrix_index = root->pmatrix_index;
	  uroot->length = uroot->back->length;
	  uroot->prob = uroot->back->prob;

	  if (!root->left)
	  {
		uroot->next = NULL;
		return uroot;
	  }

	  uroot->next = (void *)calloc(1,sizeof(pll_unetwork_node_t));
	  if (!uroot->next)
	  {
		free(uroot);
		pll_errno = PLL_ERROR_MEM_ALLOC;
		snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
		return NULL;
	  }

	  uroot->next->next = (void *)calloc(1,sizeof(pll_unetwork_node_t));
	  if (!uroot->next->next)
	  {
		free(uroot->next);
		free(uroot);
		pll_errno = PLL_ERROR_MEM_ALLOC;
		snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
		return NULL;
	  }

	  uroot->next->next->next = uroot;

	  if (root->left->is_reticulation) {
		if (root->left->first_parent == root) {
	      uroot->next->length = root->left->first_parent_length;
		  uroot->next->prob = root->left->prob;
		} else {
		  uroot->next->length = root->left->second_parent_length;
		  uroot->next->prob = 1.0 - root->left->prob;
		}
	  } else {
		uroot->next->length = root->left->length;
		uroot->next->prob = root->left->prob;
	  }
	  uroot->next->active = 1;
	  uroot->next->incoming = 0;
	  uroot->next->back = rnetwork_unroot(root->left, uroot->next, reticulation_nodes);
	  uroot->next->back->active = 1;
	  uroot->next->back->incoming = 1;

	  if (root->right->is_reticulation) {
	    if (root->right->first_parent == root) {
	  	  uroot->next->next->length = root->right->first_parent_length;
	  	  uroot->next->next->prob = root->right->prob;
	    } else {
	  	  uroot->next->next->length = root->right->second_parent_length;
	  	  uroot->next->next->prob = 1.0 - root->right->prob;
	    }
	  } else {
	  	uroot->next->next->length = root->right->length;
	  	uroot->next->next->prob = root->right->prob;
	  }
	  uroot->next->next->active = 1;
	  uroot->next->next->incoming = 0;
	  uroot->next->next->back = rnetwork_unroot(root->right, uroot->next->next, reticulation_nodes);
	  uroot->next->next->back->active = 1;
	  uroot->next->next->back->incoming = 1;
  } else { // now, we have a reticulation node
	// first, check if we already created a node for this reticulation
	if (!reticulation_nodes[root->reticulation_index]) {
      uroot = (void *)calloc(1,sizeof(pll_unetwork_node_t));
	  if (!uroot)
	  {
	    pll_errno = PLL_ERROR_MEM_ALLOC;
		snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
		return NULL;
	  }
	  uroot->back = back;
	  uroot->label = (root->label) ? xstrdup(root->label) : NULL;
	  uroot->reticulation_name = (root->reticulation_name) ? xstrdup(root->reticulation_name) : NULL;
	  uroot->reticulation_index = root->reticulation_index;
	  uroot->length = uroot->back->length;
	  uroot->prob = uroot->back->prob;

	  if (!root->child)
	  {
		free(uroot);
		pll_errno = PLL_ERROR_PARAM_INVALID;
		snprintf(pll_errmsg, 200, "The reticulation node has no child.");
		return NULL;
	  }

	  uroot->next = (void *)calloc(1,sizeof(pll_unetwork_node_t));
	  if (!uroot->next)
	  {
        free(uroot);
		pll_errno = PLL_ERROR_MEM_ALLOC;
		snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
		return NULL;
	  }

	  uroot->next->next = (void *)calloc(1,sizeof(pll_unetwork_node_t));
	  if (!uroot->next->next)
	  {
	    free(uroot->next);
		free(uroot);
		pll_errno = PLL_ERROR_MEM_ALLOC;
		snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
		return NULL;
	  }

	  reticulation_nodes[root->reticulation_index] = uroot->next->next;

	  uroot->next->next->next = uroot;

	  if (root->child->is_reticulation) {
	    if (root->child->first_parent == root) {
	      uroot->next->length = root->child->first_parent_length;
		  uroot->next->prob = root->child->prob;
	    } else {
		  uroot->next->length = root->child->second_parent_length;
		  uroot->next->prob = 1.0 - root->child->prob;
		}
	  } else {
	    uroot->next->length = root->child->length;
		uroot->next->prob = root->child->prob;
	  }
	  uroot->next->active = 1;
      uroot->next->incoming = 0;
	  uroot->next->back = rnetwork_unroot(root->child, uroot->next, reticulation_nodes);
	  uroot->next->back->active = 1;
	  uroot->next->back->incoming = 1;
	} else {
	  uroot = reticulation_nodes[root->reticulation_index];
	  uroot->back = back;
	  uroot->label = (root->label) ? xstrdup(root->label) : NULL;
	  uroot->reticulation_name = (root->reticulation_name) ? xstrdup(root->reticulation_name) : NULL;
	  uroot->length = uroot->back->length;
	  uroot->prob = uroot->back->prob;
	}
  }

  return uroot;
}

PLL_EXPORT pll_unetwork_t * pll_rnetwork_unroot(pll_rnetwork_t * network) {
  // for each node in the rnetwork, we need to create three nodes in the unetwork
  //new_network->nodes = (pll_unetwork_node_t**)malloc((network->tip_count + network->inner_tree_count * 2 + network->reticulation_count * 2) * sizeof(pll_unetwork_node_t*));
  pll_unetwork_node_t** reticulation_nodes = (pll_unetwork_node_t**)malloc(network->reticulation_count * sizeof(pll_unetwork_node_t*)); // we only need one representative per reticulation node

  /*pll_unetwork_node_t* new_root = (pll_unetwork_node_t*)malloc(sizeof(pll_unetwork_node_t*));
  new_network->vroot = new_root;
  new_root->active = 1;
  new_root->clv_index = root->clv_index;
  new_root->label = (root->label) ? xstrdup(root->label) : NULL;
  new_root->reticulation_name = (root->reticulation_name) ? xstrdup(root->reticulation_name) : NULL;*/

  pll_rnetwork_node_t * root = network->root;

  if (!root->left->left && !root->right->left)
  {
	pll_errno = PLL_ERROR_TREE_CONVERSION;
	snprintf(pll_errmsg,
			 200,
			 "Network requires at least three tips to be converted to unrooted");
	return NULL;
  }

  pll_rnetwork_node_t * new_root;

  pll_unetwork_node_t * uroot = (void *)calloc(1,sizeof(pll_unetwork_node_t));
  if (!uroot)
  {
	pll_errno = PLL_ERROR_MEM_ALLOC;
	snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
	return NULL;
  }

  uroot->next = (void *)calloc(1,sizeof(pll_unetwork_node_t));
  if (!uroot->next)
  {
	free(uroot);
	pll_errno = PLL_ERROR_MEM_ALLOC;
	snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
	return NULL;
  }

  uroot->next->next = (void *)calloc(1,sizeof(pll_unetwork_node_t));
  if (!uroot->next->next)
  {
	free(uroot->next);
	free(uroot);
	pll_errno = PLL_ERROR_MEM_ALLOC;
	snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
	return NULL;
  }

  uroot->next->next->next = uroot;

  uroot->length = root->left->length + root->right->length;
  uroot->prob = 1.0;

  /* get the first root child that has descendants and make  it the new root */
  if (root->left->left)
  {
	new_root = root->left;
	uroot->back = rnetwork_unroot(root->right,uroot, reticulation_nodes);
	/* TODO: Need to clean uroot in case of error */
	if (!uroot->back) return NULL;
  }
  else
  {
	new_root = root->right;
	uroot->back = rnetwork_unroot(root->left,uroot, reticulation_nodes);
	/* TODO: Need to clean uroot in case of error*/
	if (!uroot->back) return NULL;
  }

  uroot->label = (new_root->label) ? xstrdup(new_root->label) : NULL;
  uroot->reticulation_name = (new_root->reticulation_name) ? xstrdup(new_root->reticulation_name) : NULL;
  uroot->active = 1;
  uroot->incoming = 0;
  uroot->back->active = 1;
  uroot->back->incoming = 1;
  uroot->reticulation_index = -1;
  uroot->clv_index = root->clv_index;
  uroot->node_index = root->idx;
  uroot->scaler_index = root->scaler_index;
  uroot->pmatrix_index = root->pmatrix_index;

  uroot->next->active = 1;
  uroot->next->incoming = 0;
  uroot->next->label = uroot->label;
  uroot->next->reticulation_name = uroot->reticulation_name;
  uroot->next->reticulation_index = -1;
  uroot->next->length = new_root->left->length;
  uroot->next->prob = new_root->left->prob;
  uroot->next->back = rnetwork_unroot(new_root->left, uroot->next, reticulation_nodes);
  /* TODO: Need to clean uroot in case of error*/
  if (!uroot->next->back) return NULL;
  uroot->next->back->active = 1;
  uroot->next->back->incoming = 1;

  uroot->next->next->active = 1;
  uroot->next->next->incoming = 0;
  uroot->next->next->label = uroot->label;
  uroot->next->next->reticulation_name = uroot->reticulation_name;
  uroot->next->next->reticulation_index = -1;
  uroot->next->next->length = new_root->right->length;
  uroot->next->next->prob = new_root->right->prob;
  uroot->next->next->back = rnetwork_unroot(new_root->right, uroot->next->next, reticulation_nodes);
  /* TODO: Need to clean uroot in case of error*/
  if (!uroot->next->next->back) return NULL;
  uroot->next->next->back->active = 1;
  uroot->next->next->back->incoming = 1;

  free(reticulation_nodes);
  return pll_unetwork_wrapnetwork(uroot,0);
}

