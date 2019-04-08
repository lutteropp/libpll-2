/*
    Copyright (C) 2015 Tomas Flouri

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

static char * rnetwork_export_newick_recursive(const pll_rnetwork_node_t * root, const pll_rnetwork_node_t * taken_parent,
                                  char * (*cb_serialize)(const pll_rnetwork_node_t *))
{
  char * newick;
  int size_alloced;
  assert(root != NULL);

  if (root->is_reticulation)
  {
    if (root->first_parent->idx != taken_parent->idx)
    { // treat it as if it were a leaf
      if (cb_serialize)
      {
        newick = cb_serialize(root);
        size_alloced = strlen(newick);
      }
      else
      {
    	size_alloced = asprintf(&newick, "%s#%s:%f:%f:%f", root->label ? root->label : "", root->reticulation_name, root->support, root->length, 1.0 - root->prob);
      }
    }
    else
    { // the full subtree action
      char * subtree = rnetwork_export_newick_recursive(root->child, root, cb_serialize);
      if (subtree == NULL)
      {
    	return NULL;
      }
      if (cb_serialize)
      {
    	char * temp = cb_serialize(root);
    	size_alloced = asprintf(&newick, "(%s)%s", subtree, temp);
    	free(temp);
      }
      else
      {
    	size_alloced = asprintf(&newick, "(%s)%s#%s:%f:%f:%f", subtree, root->label ? root->label : "", root->reticulation_name, root->support, root->length, root->prob);
      }
      free(subtree);
    }
    if (size_alloced < 0)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "memory allocation during newick export failed.");
      return NULL;
    }
  }
  else
  {
    if (!(root->left) || !(root->right))
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
    else
    {
      char * subtree1 = rnetwork_export_newick_recursive(root->left, root, cb_serialize);
      if (subtree1 == NULL)
      {
        return NULL;
      }
      char * subtree2 = rnetwork_export_newick_recursive(root->right, root, cb_serialize);
      if (subtree2 == NULL)
      {
        free(subtree1);
        return NULL;
      }

      if (cb_serialize)
      {
        char * temp = cb_serialize(root);
        size_alloced = asprintf(&newick,
                                "(%s,%s)%s",
                                subtree1,
                                subtree2,
                                temp);
        free(temp);
      }
      else
      {
        size_alloced = asprintf(&newick,
                                "(%s,%s)%s:%f",
                                subtree1,
                                subtree2,
                                root->label ? root->label : "",
                                root->length);
      }
      free(subtree1);
      free(subtree2);
    }
    if (size_alloced < 0)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "memory allocation during newick export failed.");
      return NULL;
    }
  }

  return newick;
}

PLL_EXPORT char * pll_rnetwork_export_newick(const pll_rnetwork_node_t * root,
                                   char * (*cb_serialize)(const pll_rnetwork_node_t *))
{
  char * newick;
  int size_alloced;
  if (!root) return NULL;

  if (root->is_reticulation)
  {
    return PLL_FAILURE;
  }

  if (!(root->left) || !(root->right))
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
  else
  {
    char * subtree1 = rnetwork_export_newick_recursive(root->left, root, cb_serialize);
    if (subtree1 == NULL)
    {
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return NULL;
    }
    char * subtree2 = rnetwork_export_newick_recursive(root->right, root, cb_serialize);
    if (subtree2 == NULL)
    {
      free(subtree1);
      pll_errno = PLL_ERROR_MEM_ALLOC;
      snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
      return NULL;
    }

    if (cb_serialize)
    {
      char * temp = cb_serialize(root);
      size_alloced = asprintf(&newick,
                              "(%s,%s)%s",
                              subtree1,
                              subtree2,
                              temp);
      free(temp);
    }
    else
    {
      size_alloced = asprintf(&newick,
                              "(%s,%s)%s:%f;",
                              subtree1,
                              subtree2,
                              root->label ? root->label : "",
                              root->length);
    }
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0)
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "memory allocation during newick export failed");
    return NULL;
  }

  return newick;
}

static void rnetwork_tree_traverse_postorder(pll_rnetwork_node_t * node,
                                     int (*cbtrav)(pll_rnetwork_node_t *),
                                     unsigned int * index,
                                     pll_rnetwork_node_t ** outbuffer, uint64_t tree_number)
{
  if (!node)
  {
    return;
  }

  if (!node->is_reticulation)
  {
    if (!node->left && !node->right) // leaf node
	{
	  if (cbtrav(node))
	  {
	    outbuffer[*index] = node;
	    *index = *index + 1;
	  }
	  return;
	}
	if (!cbtrav(node))
	  return;


	if (node->left)
	{
      if (!node->left->is_reticulation)
      {
        rnetwork_tree_traverse_postorder(node->left, cbtrav, index, outbuffer, tree_number);
      }
      else
      {
    	int take_first_parent = (tree_number >> node->left->reticulation_index) & 1;
    	if ((node->left->first_parent->idx == node->idx && take_first_parent) || (node->left->first_parent->idx != node->idx && !take_first_parent))
    	{
    	  rnetwork_tree_traverse_postorder(node->left, cbtrav, index, outbuffer, tree_number);
    	}
      }
	}

	if (node->right)
	{
	  if (!node->right->is_reticulation)
	  {
        rnetwork_tree_traverse_postorder(node->right, cbtrav, index, outbuffer, tree_number);
	  }
	  else
	  {
		int take_first_parent = (tree_number >> node->right->reticulation_index) & 1;
		if ((node->right->first_parent->idx == node->idx && take_first_parent) || (node->right->first_parent->idx != node->idx && !take_first_parent))
		{
		  rnetwork_tree_traverse_postorder(node->right, cbtrav, index, outbuffer, tree_number);
		}
      }
	}
	outbuffer[*index] = node;
	*index = *index + 1;
  }
  else
  { // we won't add the reticulation node to the outbuffer
    rnetwork_tree_traverse_postorder(node->child, cbtrav, index, outbuffer, tree_number);
  }
}

static void rnetwork_tree_traverse_preorder(pll_rnetwork_node_t * node,
                                    int (*cbtrav)(pll_rnetwork_node_t *),
                                    unsigned int * index,
                                    pll_rnetwork_node_t ** outbuffer, uint64_t tree_number)
{
  if (!node)
  {
	return;
  }
  if (!node->is_reticulation && !node->left && !node->right)
  {
	if (cbtrav(node))
	{
	  outbuffer[*index] = node;
	  *index = *index + 1;
	}
	return;
  }

  if (!node->is_reticulation)
  {
    if (!cbtrav(node))
      return;
    outbuffer[*index] = node;
    *index = *index + 1;

	if (node->left)
	{
	  if (!node->left->is_reticulation)
	  {
	    rnetwork_tree_traverse_preorder(node->left, cbtrav, index, outbuffer, tree_number);
	  }
	  else
	  {
	    int take_first_parent = (tree_number >> node->left->reticulation_index) & 1;
	    if ((node->left->first_parent->idx == node->idx && take_first_parent) || (node->left->first_parent->idx != node->idx && !take_first_parent))
	    {
	      rnetwork_tree_traverse_preorder(node->left, cbtrav, index, outbuffer, tree_number);
	    }
	  }
	}

	if (node->right)
	{
	  if (!node->right->is_reticulation)
      {
	    rnetwork_tree_traverse_preorder(node->right, cbtrav, index, outbuffer, tree_number);
	  }
	  else
	  {
		int take_first_parent = (tree_number >> node->right->reticulation_index) & 1;
		if ((node->right->first_parent->idx == node->idx && take_first_parent) || (node->right->first_parent->idx != node->idx && !take_first_parent))
		{
		  rnetwork_tree_traverse_postorder(node->right, cbtrav, index, outbuffer, tree_number);
		}
	  }
	}
  }
  else
  { // we won't add the reticulation node to the outbuffer
	rnetwork_tree_traverse_preorder(node->child, cbtrav, index, outbuffer, tree_number);
  }
}

PLL_EXPORT int pll_rnetwork_tree_traverse(pll_rnetwork_node_t * root,
                                  int traversal,
                                  int (*cbtrav)(pll_rnetwork_node_t *),
                                  pll_rnetwork_node_t ** outbuffer,
                                  unsigned int * trav_size, uint64_t tree_number, unsigned int num_reticulations)
{
  if (tree_number >= ((unsigned int) 2 << num_reticulations))
    return PLL_FAILURE;

  *trav_size = 0;
  if (!root->left) return PLL_FAILURE;

  /* we will traverse an rooted network in the following way

           root
            /\
           /  \
        left   right

     at each node the callback function is called to decide whether we
     are going to traversing the subtree rooted at the specific node */

  if (traversal == PLL_TREE_TRAVERSE_POSTORDER)
    rnetwork_tree_traverse_postorder(root, cbtrav, trav_size, outbuffer, tree_number);
  else if (traversal == PLL_TREE_TRAVERSE_PREORDER)
    rnetwork_tree_traverse_preorder(root, cbtrav, trav_size, outbuffer, tree_number);
  else
  {
    snprintf(pll_errmsg, 200, "Invalid traversal value.");
    pll_errno = PLL_ERROR_PARAM_INVALID;
    return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}

static void rnetwork_traverse_preorder(pll_rnetwork_node_t * node,
                                    int (*cbtrav)(pll_rnetwork_node_t *),
                                    unsigned int * index,
                                    pll_rnetwork_node_t ** outbuffer)
{
  // TODO: Implement me, see https://eli.thegreenplace.net/2015/directed-graph-traversal-orderings-and-applications-to-data-flow-analysis/
}

static void rnetwork_traverse_postorder(pll_rnetwork_node_t * node,
                                    int (*cbtrav)(pll_rnetwork_node_t *),
                                    unsigned int * index,
                                    pll_rnetwork_node_t ** outbuffer)
{
  // TODO: Implement me, see https://eli.thegreenplace.net/2015/directed-graph-traversal-orderings-and-applications-to-data-flow-analysis/
}

static void rnetwork_traverse_topological(pll_rnetwork_node_t * node,
                                    int (*cbtrav)(pll_rnetwork_node_t *),
                                    unsigned int * index,
                                    pll_rnetwork_node_t ** outbuffer)
{
  // TODO: Implement me, see https://eli.thegreenplace.net/2015/directed-graph-traversal-orderings-and-applications-to-data-flow-analysis/
}

PLL_EXPORT int pll_rnetwork_traverse(pll_rnetwork_node_t * root,
                                  int traversal,
                                  int (*cbtrav)(pll_rnetwork_node_t *),
                                  pll_rnetwork_node_t ** outbuffer,
                                  unsigned int * trav_size)
{

  *trav_size = 0;
  if (!root->left) return PLL_FAILURE;

  /* we will traverse an rooted network in the following way

           root
            /\
           /  \
        left   right

     at each node the callback function is called to decide whether we
     are going to traversing the subtree rooted at the specific node */

  if (traversal == PLL_NETWORK_TRAVERSE_POSTORDER)
    rnetwork_traverse_postorder(root, cbtrav, trav_size, outbuffer);
  else if (traversal == PLL_NETWORK_TRAVERSE_PREORDER)
    rnetwork_traverse_preorder(root, cbtrav, trav_size, outbuffer);
  else if (traversal == PLL_NETWORK_TRAVERSE_TOPOLOGICAL)
	rnetwork_traverse_topological(root, cbtrav, trav_size, outbuffer);
  else
  {
    snprintf(pll_errmsg, 200, "Invalid traversal value.");
    pll_errno = PLL_ERROR_PARAM_INVALID;
    return PLL_FAILURE;
  }

  return PLL_SUCCESS;
}

PLL_EXPORT pll_rtree_t * pll_rnetwork_extract_rtree(pll_rnetwork_node_t * root, uint64_t tree_number, unsigned int num_reticulations)
{
  // TODO: Implement me.
	return NULL;
}

