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
	assert(node);
	unsigned int cnt = 0;
	if (!node->incoming) {
		cnt++;
	}
	pll_unetwork_node_t * snode = node->next;
	while (snode && snode != node) {
		if (!snode->incoming) {
			cnt++;
		}
		snode = snode->next;
	}
	return cnt;
}

PLL_EXPORT unsigned int pll_unetwork_count_active_outgoing(const pll_unetwork_node_t * node) {
	assert(node);
	unsigned int cnt = 0;
	if (!node->incoming && node->active) {
		cnt++;
	}
	pll_unetwork_node_t * snode = node->next;
	while (snode && snode != node) {
		if (!snode->incoming && snode->active) {
			cnt++;
		}
		snode = snode->next;
	}
	return cnt;
}

unsigned int count_incoming(const pll_unetwork_node_t * node) {
	assert(node);
	unsigned int cnt = 0;
	if (node->incoming) {
		cnt++;
	}
	pll_unetwork_node_t * snode = node->next;
	while (snode && snode != node) {
		if (snode->incoming) {
			cnt++;
		}
		snode = snode->next;
	}
	return cnt;
}

PLL_EXPORT unsigned int pll_unetwork_count_active_incoming(const pll_unetwork_node_t * node) {
	assert(node);
	unsigned int cnt = 0;
	if (node->incoming && node->active) {
		cnt++;
	}
	pll_unetwork_node_t * snode = node->next;
	while (snode && snode != node) {
		if (snode->incoming && snode->active) {
			cnt++;
		}
		snode = snode->next;
	}
	return cnt;
}

PLL_EXPORT int pll_unetwork_is_inner_tree(const pll_unetwork_node_t * node) {
	assert(node);
	unsigned int cnt_out = count_outgoing(node);
	return (cnt_out > 1);
}

PLL_EXPORT int pll_unetwork_is_reticulation(const pll_unetwork_node_t * node) {
	assert(node);
	unsigned int cnt_in = count_incoming(node);
	return (cnt_in > 1);
}

PLL_EXPORT int pll_unetwork_is_leaf(const pll_unetwork_node_t * node) {
	assert(node);
	unsigned int cnt_out = count_outgoing(node);
	return (cnt_out == 0);
}

PLL_EXPORT int pll_unetwork_is_root(const pll_unetwork_node_t * node) {
	assert(node);
	unsigned int cnt_in = count_incoming(node);
	return (cnt_in == 0);
}

PLL_EXPORT pll_unetwork_node_t * pll_unetwork_get_reticulation_child(const pll_unetwork_node_t * node) {
	assert(pll_unetwork_is_reticulation(node));
	while (node && node->incoming) {
		node = node->next;
	}
	assert(node && !node->incoming && node->active);
	return node->back;
}

PLL_EXPORT pll_unetwork_node_t * pll_unetwork_get_active_parent(const pll_unetwork_node_t* node) {
	assert(pll_unetwork_count_active_incoming(node) == 1);
	while (node && (!node->incoming || !node->active)) {
		node = node->next;
	}
	assert(node && node->incoming && node->active);
	return node->back;
}

PLL_EXPORT int pll_unetwork_get_tree_children(const pll_unetwork_node_t * node, pll_unetwork_node_t ** child1, pll_unetwork_node_t ** child2) {
	assert(node);
	assert(!pll_unetwork_is_reticulation(node));
	pll_unetwork_node_t * snode = node;
	do {
		if (!snode->incoming) {
			if (!*child1) {
				*child1 = snode->back;
			} else if (!*child2) {
				*child2 = snode->back;
			} else {
				return PLL_FAILURE;
			}
		}
		snode = snode->next;
	} while (snode && snode != node);
	return PLL_SUCCESS;
}

static char * newick_unetwork_recurse(const pll_unetwork_node_t * root,
                                    char * (*cb_serialize)(const pll_unetwork_node_t *),
                                    int level, int * visited_reticulations)
{
  char * newick;
  int size_alloced = 0;
  assert(root != NULL);

  if (pll_unetwork_is_reticulation(root) && visited_reticulations[root->reticulation_index]) { // we already encountered this node, act like if it is a dead end.
	if (cb_serialize)
	{
	  // TODO: does this work?
	  newick = cb_serialize(root);
	  size_alloced = strlen(newick);
	}
	else
	{
	  if (root->reticulation_name)
	  {
	    size_alloced = asprintf(&newick, "%s#%s:%f:%f:%f", root->label ? root->label : "", root->reticulation_name, root->length, root->support, root->prob);
	  }
	  else
	  {
		size_alloced = asprintf(&newick, "%s#%d:%f:%f:%f", root->label ? root->label : "", root->reticulation_index, root->length, root->support, root->prob);
	  }
	}
  }
  else if (!root->next) // leaf node
  {
    if (cb_serialize)
    {
      newick = cb_serialize(root);
      size_alloced = strlen(newick);
    }
    else
    {
      size_alloced = asprintf(&newick, "%s:%f", root->label, root->length);
    }
  }
  else // inner node
  {
    const pll_unetwork_node_t * start = root->next;
    const pll_unetwork_node_t * snode = start;
    char * cur_newick;
    do
    {
      if (!snode->incoming) {
        char * subnetwork = newick_unetwork_recurse(snode->back, cb_serialize, level+1, visited_reticulations);
        if (subnetwork == NULL)
        {
          pll_errno = PLL_ERROR_MEM_ALLOC;
          snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
          return NULL;
        }

        if (snode == start)
        {
          cur_newick = subnetwork;
        }
        else
        {
          char * temp = cur_newick;
          size_alloced = asprintf(&cur_newick,
                                  "%s,%s",
                                  temp,
                                  subnetwork);
          free(temp);
          free(subnetwork);
        }
      }
      snode = snode->next;
    }
    while(snode != root);

    if (level > 0)
    {
      if (pll_unetwork_is_reticulation(root))
      {
    	if (cb_serialize)
		{
		  char * temp = cb_serialize(root);
		  size_alloced = asprintf(&newick,
		  						"(%s)%s",
		  						cur_newick,
		  						temp);
		  free(temp);
		}
		else
		{
		  if (root->reticulation_name)
		  {
		    size_alloced = asprintf(&newick,
			  					  "(%s)%s#%s:%f:%f:%f",
				  				  cur_newick,
					  			  root->label ? root->label : "",
						  		  root->reticulation_name,
							  	  root->length,
								  root->support,
								  root->prob);
		  }
		  else
		  {
		    size_alloced = asprintf(&newick,
								  "(%s)%s#%d:%f:%f:%f",
								  cur_newick,
								  root->label ? root->label : "",
								  root->reticulation_index,
								  root->length,
								  root->support,
					              root->prob);
		  }
		}
		free(cur_newick);
        visited_reticulations[root->reticulation_index] = 1;
      }
      else
      {
        if (cb_serialize)
        {
          char * temp = cb_serialize(root);
          size_alloced = asprintf(&newick,
                                  "(%s)%s",
                                  cur_newick,
                                  temp);
          free(temp);
        }
        else
        {
          size_alloced = asprintf(&newick,
                                  "(%s)%s:%f",
                                  cur_newick,
                                  root->label ? root->label : "",
                                  root->length);
        }
        free(cur_newick);
      }
    }
    else
      newick = cur_newick;
  }

  if (size_alloced < 0)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "memory allocation during newick export failed");
    return NULL;
  }

  return newick;
}

char * unetwork_export_newick(const pll_unetwork_t * network,
                           int export_rooted,
                           double root_brlen,
                           char * (*cb_serialize)(const pll_unetwork_node_t *))
{
  char * newick;
  char * subtree1;
  char * subtree2;
  int size_alloced;

  pll_unetwork_node_t * root = network->vroot;

  if (!root) return NULL;

  int* visited_reticulations = (int*) calloc(network->reticulation_count, sizeof(int));

  if (!root->next) root = root->back;

  if (export_rooted)
  {
    assert(!cb_serialize);

    subtree1 = newick_unetwork_recurse(root->back, cb_serialize, 1, visited_reticulations);
    subtree2 = newick_unetwork_recurse(root, cb_serialize, 0, visited_reticulations);

    size_alloced = asprintf(&newick,
                            "(%s,(%s)%s:%f):0.0;",
                            subtree1,
                            subtree2,
                            root->label ? root->label : "",
                            root_brlen);
  }
  else
  {
    subtree1 = newick_unetwork_recurse(root->back, cb_serialize, 1, visited_reticulations);
    subtree2 = newick_unetwork_recurse(root, cb_serialize, 0, visited_reticulations);

    size_alloced = asprintf(&newick,
                            "(%s,%s)%s:0.0;",
                            subtree1,
                            subtree2,
                            root->label ? root->label : "");
  }

  free(subtree1);
  free(subtree2);

  free(visited_reticulations);

  if (size_alloced < 0)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "memory allocation during newick export failed");
    return NULL;
  }

//  printf("newick: %s\n", newick);

  return (newick);
}

PLL_EXPORT char * pll_unetwork_export_newick(const pll_unetwork_t * network,
                                   char * (*cb_serialize)(const pll_unetwork_node_t *)) {
	if (network->reticulation_count == 0) {
	  return unetwork_export_newick(network, 0, 0, cb_serialize);
	} else {
	  return unetwork_export_newick(network, 1, 0, cb_serialize);
	}
}

static void unetwork_traverse_recursive(pll_unetwork_node_t * node,
		                                int traversal,
		                                int (*cbtrav)(pll_unetwork_node_t *),
		                                unsigned int * index,
								        pll_unetwork_node_t ** outbuffer,
	                                    int* visited_reticulations) // TODO: Is this correct?
{
  if (!node || (pll_unetwork_is_reticulation(node) && visited_reticulations[node->reticulation_index]) || !cbtrav(node)) {
	  return;
  }

  // if preorder traversal, first deal with the node.
  if (traversal == PLL_NETWORK_TRAVERSE_PREORDER) {
	assert(node->incoming || pll_unetwork_is_root(node));
	outbuffer[*index] = node;
	*index = *index + 1;
	if (pll_unetwork_is_reticulation(node))
	{
	  visited_reticulations[node->reticulation_index] = 1;
	}
  }

  pll_unetwork_node_t * snode = node;
  do
  {
	if (!snode->incoming && snode->active) // outgoing edge
	{
		unetwork_traverse_recursive(snode->back, traversal, cbtrav, index, outbuffer, visited_reticulations);
	}
	snode = snode->next;
  } while (snode && snode != node);

  // if postorder traversal, first deal with the children
  if (traversal == PLL_NETWORK_TRAVERSE_POSTORDER) {
	assert(node->incoming || pll_unetwork_is_root(node));
	outbuffer[*index] = node;
	*index = *index + 1;
	if (pll_unetwork_is_reticulation(node))
	{
	  visited_reticulations[node->reticulation_index] = 1;
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
    pll_unetwork_node_t * snode = network->reticulation_nodes[i];
    pll_unetwork_node_t * first_incoming = NULL;
    pll_unetwork_node_t * second_incoming = NULL;
    do {
      if (snode->incoming) {
    	if (!first_incoming) {
    	  first_incoming = snode;
    	} else {
    	  second_incoming = snode;
    	}
      }
      snode = snode->next;
    } while (snode && snode != network->reticulation_nodes[i]);
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

PLL_EXPORT int pll_unetwork_traverse(pll_unetwork_node_t * root,
                                  int traversal,
                                  int (*cbtrav)(pll_unetwork_node_t *),
                                  pll_unetwork_node_t ** outbuffer,
                                  unsigned int * trav_size) {
  *trav_size = 0;

  if (traversal == PLL_NETWORK_TRAVERSE_POSTORDER ||
	  traversal == PLL_NETWORK_TRAVERSE_PREORDER)
  {
	int* visited_reticulations = (int*) calloc(MAX_RETICULATION_COUNT, sizeof(int));
	unetwork_traverse_recursive(root, traversal, cbtrav, trav_size, outbuffer, visited_reticulations);
	free(visited_reticulations);
  }
  else
  {
	snprintf(pll_errmsg, 200, "Invalid traversal value.");
	pll_errno = PLL_ERROR_PARAM_INVALID;
	return PLL_FAILURE;
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
  return pll_unetwork_traverse(network->vroot, traversal, cbtrav, outbuffer, trav_size);
}

PLL_EXPORT void pll_unetwork_create_operations(pll_unetwork_node_t * const* trav_buffer,
                                            unsigned int trav_buffer_size,
                                            double * branches,
                                            unsigned int * pmatrix_indices,
                                            pll_operation_t * ops,
                                            unsigned int * matrix_count,
                                            unsigned int * ops_count) {
  const pll_unetwork_node_t * node;
  unsigned int i;

  *ops_count = 0;
  if (matrix_count)
	*matrix_count = 0;

  for (i = 0; i < trav_buffer_size; ++i)
  {
	node = trav_buffer[i];

	/* if the current node is the second end-point of the edge
	shared with the root node, then do not add the edge to the
	list as it will be added in the end (avoid duplicate edges
	in the list) */
	if (node != trav_buffer[trav_buffer_size - 1]->back)
	{
	  if (branches)
		*branches++ = node->length;
	  if (pmatrix_indices)
		*pmatrix_indices++ = node->pmatrix_index;
	  if (matrix_count)
		*matrix_count = *matrix_count + 1;
	}

	if (node->next)
	{
	  ops[*ops_count].parent_clv_index = node->clv_index;
	  ops[*ops_count].parent_scaler_index = node->scaler_index;

	  ops[*ops_count].child1_clv_index = node->next->back->clv_index;
	  ops[*ops_count].child1_scaler_index = node->next->back->scaler_index;
	  ops[*ops_count].child1_matrix_index = node->next->back->pmatrix_index;

	  ops[*ops_count].child2_clv_index = node->next->next->back->clv_index;
	  ops[*ops_count].child2_scaler_index = node->next->next->back->scaler_index;
	  ops[*ops_count].child2_matrix_index = node->next->next->back->pmatrix_index;

	  *ops_count = *ops_count + 1;
	}
  }
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

  if (pll_unetwork_count_active_incoming(node) > 1) {
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

      if (network->binary && subnodes > 3 && !pll_unetwork_is_reticulation(snode))
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

      if (network->binary && subnodes > 3 && !pll_unetwork_is_reticulation(snode))
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

PLL_EXPORT pll_unetwork_t * pll_unetwork_clone(const pll_unetwork_t * network) {
   /* choose the last inner node as the starting point of the clone. It does not
	 really matter which node to choose, but since the newick parser places the
	 root node at the end of the list, we use the same notation here */

  pll_unetwork_t * cloned_network = (pll_unetwork_t *) malloc(sizeof(pll_unetwork_t));
  cloned_network->binary = network->binary;
  cloned_network->edge_count = network->edge_count;
  cloned_network->inner_tree_count = network->inner_tree_count;
  cloned_network->reticulation_count = network->reticulation_count;
  cloned_network->tip_count = network->tip_count;
  cloned_network->tree_edge_count = network->tree_edge_count;
  unsigned int total_node_count = network->inner_tree_count + network->reticulation_count + network->tip_count;
  cloned_network->nodes = (pll_unetwork_node_t **) malloc(total_node_count * sizeof(pll_unetwork_node_t *));
  cloned_network->reticulation_nodes = (pll_unetwork_node_t **) malloc(network->reticulation_count * sizeof(pll_unetwork_node_t *));

  unsigned int total_mini_nodes_count = network->tip_count + 3 * (network->reticulation_count + network->inner_tree_count);
  pll_unetwork_node_t ** orig_node_mappings = (pll_unetwork_node_t **) malloc(total_mini_nodes_count * sizeof(pll_unetwork_t *));
  pll_unetwork_node_t ** cloned_node_mappings = (pll_unetwork_node_t **) malloc(total_mini_nodes_count * sizeof(pll_unetwork_t *));
  unsigned int i;
  for (i = 0; i < total_node_count; ++i) {
	  // copy each node, as well as its next ones
	  pll_unetwork_node_t * node = network->nodes[i];

	  // first, the "original" node
	  pll_unetwork_node_t * new_node = (pll_unetwork_node_t *) malloc(sizeof(pll_unetwork_node_t));
	  if (!new_node) {
		pll_errno = PLL_ERROR_MEM_ALLOC;
		snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
		return NULL;
	  }
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

	  orig_node_mappings[node->node_index] = node;
	  cloned_node_mappings[node->node_index] = new_node;

	  // now, the other links
	  pll_unetwork_node_t * snode = node->next;
	  while (snode && snode != node) {
          pll_unetwork_node_t * new_snode = (pll_unetwork_node_t *) malloc(sizeof(pll_unetwork_node_t));
          if (!new_snode) {
		    pll_errno = PLL_ERROR_MEM_ALLOC;
		    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
		    return NULL;
		  }
		  memcpy(new_snode, snode, sizeof(pll_unetwork_node_t));
		  new_snode->label = new_node->label;
		  new_snode->reticulation_name = new_node->reticulation_name;
          orig_node_mappings[snode->node_index] = snode;
          cloned_node_mappings[snode->node_index] = new_snode;

		  snode = snode->next;
	  }
  }
  // now, deal with the connections
  for (i = 0; i < total_mini_nodes_count; ++i) {
	  assert(cloned_node_mappings[i]);
	  assert(orig_node_mappings[i]);
	  cloned_node_mappings[i]->back = cloned_node_mappings[orig_node_mappings[i]->back->node_index];
	  if (orig_node_mappings[i]->next) { // non-leaf node
	    cloned_node_mappings[i]->next = cloned_node_mappings[orig_node_mappings[i]->next->node_index];
	  }
  }
  // now, deal with the reticulations
  for (i = 0; i < network->reticulation_count; ++i) {
	  cloned_network->reticulation_nodes[i] = cloned_node_mappings[network->reticulation_nodes[i]->node_index];
  }
  // now, deal with the nodes
  for (i = 0; i < network->inner_tree_count + network->reticulation_count + network->tip_count; ++i) {
	  cloned_network->nodes[i] = cloned_node_mappings[network->nodes[i]->node_index];
  }
  // now, deal with the root
  cloned_network->vroot = cloned_node_mappings[network->vroot->node_index];
  free(cloned_node_mappings);
  free(orig_node_mappings);

 /* pll_unetwork_node_t * root = pll_unetwork_graph_clone(network->vroot);

  if (network->binary)
	return pll_unetwork_wrapnetwork(root, network->tip_count);
  else
	return pll_unetwork_wrapnetwork_multi(root, network->tip_count, network->inner_tree_count, network->reticulation_count);*/
  return cloned_network;
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
	  uroot->reticulation_name = NULL;
	  uroot->reticulation_index = -1;
	  uroot->length = uroot->back->length;
	  uroot->prob = uroot->back->prob;

	  if (!root->left && !root->right)
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
		  uroot->next->prob = root->left->first_parent_prob;
		} else {
		  uroot->next->length = root->left->second_parent_length;
		  uroot->next->prob = root->left->second_parent_prob;
		}
	  } else {
		uroot->next->length = root->left->length;
		uroot->next->prob = 1.0;
	  }
	  uroot->next->label = uroot->label;
	  uroot->next->reticulation_index = uroot->reticulation_index;
	  uroot->next->reticulation_name = uroot->reticulation_name;
	  uroot->next->active = 1;
	  uroot->next->incoming = 0;
	  uroot->next->back = rnetwork_unroot(root->left, uroot->next, reticulation_nodes);
	  uroot->next->back->active = 1;
	  uroot->next->back->incoming = 1;

	  if (root->right->is_reticulation) {
	    if (root->right->first_parent == root) {
	  	  uroot->next->next->length = root->right->first_parent_length;
	  	  uroot->next->next->prob = root->right->first_parent_prob;
	    } else {
	  	  uroot->next->next->length = root->right->second_parent_length;
	  	  uroot->next->next->prob = root->right->second_parent_prob;
	    }
	  } else {
	  	uroot->next->next->length = root->right->length;
	  	uroot->next->next->prob = 1.0;
	  }
	  uroot->next->next->label = uroot->label;
	  uroot->next->next->reticulation_index = uroot->reticulation_index;
	  uroot->next->next->reticulation_name = uroot->reticulation_name;

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
		  uroot->next->prob = root->child->first_parent_prob;
	    } else {
		  uroot->next->length = root->child->second_parent_length;
		  uroot->next->prob = 1.0 - root->child->second_parent_prob;
		}
	  } else {
	    uroot->next->length = root->child->length;
		uroot->next->prob = 1.0;
	  }
	  uroot->next->label = uroot->label;
	  uroot->next->reticulation_name = uroot->reticulation_name;
	  uroot->next->reticulation_index = uroot->reticulation_index;

	  uroot->next->active = 1;
      uroot->next->incoming = 0;
	  uroot->next->back = rnetwork_unroot(root->child, uroot->next, reticulation_nodes);
	  uroot->next->back->active = 1;
	  uroot->next->back->incoming = 1;

	  uroot->next->next->label = uroot->label;
	  uroot->next->next->reticulation_name = uroot->reticulation_name;
	  uroot->next->next->reticulation_index = uroot->reticulation_index;
	} else {
	  uroot = reticulation_nodes[root->reticulation_index];
	  uroot->back = back;
	  uroot->prob = uroot->back->prob;
	  uroot->length = uroot->back->length;
	}
  }

  return uroot;
}

PLL_EXPORT pll_unetwork_t * pll_rnetwork_unroot(pll_rnetwork_t * network) {
  // for each node in the rnetwork, we need to create three nodes in the unetwork... except for the leaves.
  pll_unetwork_node_t** reticulation_nodes = (pll_unetwork_node_t**)malloc(network->reticulation_count * sizeof(pll_unetwork_node_t*)); // we only need one representative per reticulation node

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

  /* get the first root child that has descendants and make it the new root */
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
  uroot->reticulation_name = NULL;
  uroot->active = 1;
  uroot->incoming = 0;
  uroot->back->active = 1;
  uroot->back->incoming = 1;
  uroot->reticulation_index = -1;

  uroot->next->active = 1;
  uroot->next->incoming = 0;
  uroot->next->label = uroot->label;
  uroot->next->reticulation_name = uroot->reticulation_name;
  uroot->next->reticulation_index = -1;
  if (new_root->left->is_reticulation)
  {
    if (new_root->left->first_parent == new_root)
    {
      uroot->next->prob = new_root->left->first_parent_prob;
      uroot->next->length = new_root->left->first_parent_length;
    } else
    {
      uroot->next->prob = new_root->left->second_parent_prob;
      uroot->next->length = new_root->left->second_parent_length;
    }
  }
  else
  {
  	uroot->next->length = new_root->left->length;
    uroot->next->prob = 1.0;
  }
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

  if (new_root->right->is_reticulation)
  {
    if (new_root->right->first_parent == new_root)
    {
      uroot->next->next->prob = new_root->right->first_parent_prob;
      uroot->next->next->length = new_root->right->first_parent_length;
    } else
    {
      uroot->next->next->prob = new_root->right->second_parent_prob;
      uroot->next->next->length = new_root->right->second_parent_length;
    }
  }
  else
  {
	uroot->next->next->length = new_root->right->length;
	uroot->next->next->prob = 1.0;
  }
  uroot->next->next->back = rnetwork_unroot(new_root->right, uroot->next->next, reticulation_nodes);
  /* TODO: Need to clean uroot in case of error*/
  if (!uroot->next->next->back) return NULL;
  uroot->next->next->back->active = 1;
  uroot->next->next->back->incoming = 1;

  free(reticulation_nodes);
  pll_unetwork_t * res = pll_unetwork_wrapnetwork(uroot,0);
  // re-wire pointer to fit the utree conventions...
  res->vroot = res->vroot->next;
  return res;
}

PLL_EXPORT void pll_unetwork_create_pars_buildops(pll_unetwork_node_t * const* trav_buffer,
                                               unsigned int trav_buffer_size,
                                               pll_pars_buildop_t * ops,
                                               unsigned int * ops_count)
{
  const pll_unetwork_node_t * node;
  unsigned int i;

  *ops_count = 0;

  for (i = 0; i < trav_buffer_size; ++i)
  {
    node = trav_buffer[i];

    if (node->next)
    {
      ops[*ops_count].parent_score_index = node->clv_index;
      ops[*ops_count].child1_score_index = node->next->back->clv_index;
      ops[*ops_count].child2_score_index = node->next->next->back->clv_index;

      *ops_count = *ops_count + 1;
    }
  }
}

