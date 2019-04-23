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
		char * (*cb_serialize)(const pll_rnetwork_node_t *)) {
	char * newick;
	int size_alloced;
	assert(root != NULL);

	if (root->is_reticulation) {
		if (root->first_parent->idx != taken_parent->idx) { // treat it as if it were a leaf
			if (cb_serialize) {
				newick = cb_serialize(root);
				size_alloced = strlen(newick);
			} else {
				size_alloced = asprintf(&newick, "%s#%s:%f:%f:%f", root->label ? root->label : "", root->reticulation_name, root->support,
						root->second_parent_length, 1.0 - root->prob);
			}
		} else { // the full subtree action
			char * subtree = rnetwork_export_newick_recursive(root->child, root, cb_serialize);
			if (subtree == NULL) {
				return NULL;
			}
			if (cb_serialize) {
				char * temp = cb_serialize(root);
				size_alloced = asprintf(&newick, "(%s)%s", subtree, temp);
				free(temp);
			} else {
				size_alloced = asprintf(&newick, "(%s)%s#%s:%f:%f:%f", subtree, root->label ? root->label : "", root->reticulation_name,
						root->support, root->first_parent_length, root->prob);
			}
			free(subtree);
		}
		if (size_alloced < 0) {
			pll_errno = PLL_ERROR_MEM_ALLOC;
			snprintf(pll_errmsg, 200, "memory allocation during newick export failed.");
			return NULL;
		}
	} else {
		if (!(root->left) || !(root->right)) {
			if (cb_serialize) {
				newick = cb_serialize(root);
				size_alloced = strlen(newick);
			} else {
				size_alloced = asprintf(&newick, "%s:%f", root->label, root->length);
			}
		} else {
			char * subtree1 = rnetwork_export_newick_recursive(root->left, root, cb_serialize);
			if (subtree1 == NULL) {
				return NULL;
			}
			char * subtree2 = rnetwork_export_newick_recursive(root->right, root, cb_serialize);
			if (subtree2 == NULL) {
				free(subtree1);
				return NULL;
			}

			if (cb_serialize) {
				char * temp = cb_serialize(root);
				size_alloced = asprintf(&newick, "(%s,%s)%s", subtree1, subtree2, temp);
				free(temp);
			} else {
				size_alloced = asprintf(&newick, "(%s,%s)%s:%f", subtree1, subtree2, root->label ? root->label : "", root->length);
			}
			free(subtree1);
			free(subtree2);
		}
		if (size_alloced < 0) {
			pll_errno = PLL_ERROR_MEM_ALLOC;
			snprintf(pll_errmsg, 200, "memory allocation during newick export failed.");
			return NULL;
		}
	}

	return newick;
}

PLL_EXPORT char * pll_rnetwork_export_newick(const pll_rnetwork_node_t * root, char * (*cb_serialize)(const pll_rnetwork_node_t *)) {
	char * newick;
	int size_alloced;
	if (!root)
		return NULL;

	if (root->is_reticulation) {
		return PLL_FAILURE;
	}

	if (!(root->left) || !(root->right)) {
		if (cb_serialize) {
			newick = cb_serialize(root);
			size_alloced = strlen(newick);
		} else {
			size_alloced = asprintf(&newick, "%s:%f", root->label, root->length);
		}
	} else {
		char * subtree1 = rnetwork_export_newick_recursive(root->left, root, cb_serialize);
		if (subtree1 == NULL) {
			pll_errno = PLL_ERROR_MEM_ALLOC;
			snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
			return NULL;
		}
		char * subtree2 = rnetwork_export_newick_recursive(root->right, root, cb_serialize);
		if (subtree2 == NULL) {
			free(subtree1);
			pll_errno = PLL_ERROR_MEM_ALLOC;
			snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
			return NULL;
		}

		if (cb_serialize) {
			char * temp = cb_serialize(root);
			size_alloced = asprintf(&newick, "(%s,%s)%s", subtree1, subtree2, temp);
			free(temp);
		} else {
			size_alloced = asprintf(&newick, "(%s,%s)%s:%f;", subtree1, subtree2, root->label ? root->label : "", root->length);
		}
		free(subtree1);
		free(subtree2);
	}
	if (size_alloced < 0) {
		pll_errno = PLL_ERROR_MEM_ALLOC;
		snprintf(pll_errmsg, 200, "memory allocation during newick export failed");
		return NULL;
	}

	return newick;
}

PLL_EXPORT double pll_rnetwork_reticulation_logprob(pll_rnetwork_t * network, uint64_t tree_number) {
	double res = 0;
	if (tree_number >= ((unsigned int) 2 << network->reticulation_count))
		return PLL_FAILURE;
	unsigned int i;
	// TODO: This can be made more efficient by storing pointers to the reticulation nodes
	for (i = 0; i < network->tip_count + network->inner_tree_count + network->reticulation_count; ++i) {
		if (network->nodes[i]->is_reticulation) {
			if (pll_rnetwork_can_go_tree(network->nodes[i]->first_parent, network->nodes[i], tree_number)) {
				res += log(network->nodes[i]->prob);
			} else {
				res += log(1.0 - network->nodes[i]->prob);
			}
		}
	}
	return res;
}

PLL_EXPORT double pll_rnetwork_reticulation_prob(pll_rnetwork_t * network, uint64_t tree_number) {
	double res = 1.0;
	if (tree_number >= ((unsigned int) 2 << network->reticulation_count))
		return PLL_FAILURE;
	unsigned int i;
	// TODO: This can be made more efficient by storing pointers to the reticulation nodes
	for (i = 0; i < network->tip_count + network->inner_tree_count + network->reticulation_count; ++i) {
		if (network->nodes[i]->is_reticulation) {
			if (pll_rnetwork_can_go_tree(network->nodes[i]->first_parent, network->nodes[i], tree_number)) {
				res *= network->nodes[i]->prob;
			} else {
				res *= 1.0 - network->nodes[i]->prob;
			}
		}
	}
	return res;
}

PLL_EXPORT int pll_rnetwork_can_go_tree(pll_rnetwork_node_t * parent, pll_rnetwork_node_t * child, uint64_t tree_number) {
	if (!child->is_reticulation) {
		return 1;
	}
	int take_first_parent = (tree_number >> child->reticulation_index) & 1;
	if ((child->first_parent->idx == parent->idx && take_first_parent) || (child->first_parent->idx != parent->idx && !take_first_parent)) {
		return 1;
	} else {
		return 0;
	}
}

static void rnetwork_tree_traverse_postorder(pll_rnetwork_node_t * node, int (*cbtrav)(pll_rnetwork_node_t *), unsigned int * index,
		pll_rnetwork_node_t ** outbuffer, uint64_t tree_number, int* dead) {
	if (!node) {
		return;
	}

	if (!node->is_reticulation) {
		if (!node->left && !node->right) // leaf node
				{
			if (cbtrav(node)) {
				outbuffer[*index] = node;
				*index = *index + 1;
			}
			return;
		}
		if (!cbtrav(node))
			return;

		unsigned int living_children_count = 0;

		if (node->left) {
			if (!node->left->is_reticulation) {
				rnetwork_tree_traverse_postorder(node->left, cbtrav, index, outbuffer, tree_number, dead);
				if (!dead[node->left->idx]) {
					living_children_count++;
				}
			} else {
				if (pll_rnetwork_can_go_tree(node, node->left, tree_number)) {
					rnetwork_tree_traverse_postorder(node->left, cbtrav, index, outbuffer, tree_number, dead);
					if (!dead[node->left->idx]) {
						living_children_count++;
					}
				}
			}
		}

		if (node->right) {
			if (!node->right->is_reticulation) {
				rnetwork_tree_traverse_postorder(node->right, cbtrav, index, outbuffer, tree_number, dead);
				if (!dead[node->right->idx]) {
					living_children_count++;
				}
			} else {
				if (pll_rnetwork_can_go_tree(node, node->right, tree_number)) {
					rnetwork_tree_traverse_postorder(node->right, cbtrav, index, outbuffer, tree_number, dead);
					if (!dead[node->right->idx]) {
						living_children_count++;
					}
				}
			}
		}

		if (living_children_count == 0) {
			dead[node->idx] = 1;
		} else if (living_children_count > 1) // both children are alive.
				{
			outbuffer[*index] = node;
			*index = *index + 1;
		}
	} else { // we won't add the reticulation node to the outbuffer
		rnetwork_tree_traverse_postorder(node->child, cbtrav, index, outbuffer, tree_number, dead);
	}
}

PLL_EXPORT int pll_rnetwork_tree_traverse(pll_rnetwork_t * network, int traversal, int (*cbtrav)(pll_rnetwork_node_t *),
		pll_rnetwork_node_t ** outbuffer, unsigned int * trav_size, uint64_t tree_number) {
	if (tree_number >= ((unsigned int) 2 << network->reticulation_count))
		return PLL_FAILURE;

	*trav_size = 0;
	if (!network->root->left)
		return PLL_FAILURE;

	unsigned int total_node_count = network->tip_count + network->inner_tree_count + network->reticulation_count;
	int* dead = (int*) calloc(total_node_count, sizeof(int));

	/* we will traverse an rooted network in the following way

	 root
	 /\
           /  \
        left   right

	 at each node the callback function is called to decide whether we
	 are going to traversing the subtree rooted at the specific node */

	if (traversal == PLL_TREE_TRAVERSE_POSTORDER)
		rnetwork_tree_traverse_postorder(network->root, cbtrav, trav_size, outbuffer, tree_number, dead);
	else {
		snprintf(pll_errmsg, 200, "Invalid traversal value.");
		pll_errno = PLL_ERROR_PARAM_INVALID;
		free(dead);
		return PLL_FAILURE;
	}

	free(dead);
	return PLL_SUCCESS;
}

static void rnetwork_traverse_preorder(pll_rnetwork_node_t * node, int (*cbtrav)(pll_rnetwork_node_t *), unsigned int * index,
		pll_rnetwork_node_t ** outbuffer) {
	// TODO: Implement me, see https://eli.thegreenplace.net/2015/directed-graph-traversal-orderings-and-applications-to-data-flow-analysis/
}

static void rnetwork_traverse_postorder(pll_rnetwork_node_t * node, int (*cbtrav)(pll_rnetwork_node_t *), unsigned int * index,
		pll_rnetwork_node_t ** outbuffer) {
	// TODO: Implement me, see https://eli.thegreenplace.net/2015/directed-graph-traversal-orderings-and-applications-to-data-flow-analysis/
}

static void rnetwork_traverse_topological(pll_rnetwork_node_t * node, int (*cbtrav)(pll_rnetwork_node_t *), unsigned int * index,
		pll_rnetwork_node_t ** outbuffer) {
	// TODO: Implement me, see https://eli.thegreenplace.net/2015/directed-graph-traversal-orderings-and-applications-to-data-flow-analysis/
}

PLL_EXPORT int pll_rnetwork_traverse(pll_rnetwork_node_t * root, int traversal, int (*cbtrav)(pll_rnetwork_node_t *),
		pll_rnetwork_node_t ** outbuffer, unsigned int * trav_size) {

	*trav_size = 0;
	if (!root->left)
		return PLL_FAILURE;

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
	else {
		snprintf(pll_errmsg, 200, "Invalid traversal value.");
		pll_errno = PLL_ERROR_PARAM_INVALID;
		return PLL_FAILURE;
	}

	return PLL_SUCCESS;
}

/* TODO: Memory allocation checks were not implemented in this function!!! */
static pll_rnetwork_node_t * clone_node(const pll_rnetwork_node_t * node) {
	pll_rnetwork_node_t * new_node = (pll_rnetwork_node_t *) malloc(sizeof(pll_rnetwork_node_t));
	memcpy(new_node, node, sizeof(pll_rnetwork_node_t));

	if (node->label) {
		new_node->label = (char *) malloc(strlen(node->label) + 1);
		strcpy(new_node->label, node->label);
	}
	if (node->reticulation_name) {
		new_node->reticulation_name = (char *) malloc(strlen(node->reticulation_name) + 1);
		strcpy(new_node->reticulation_name, node->reticulation_name);
	}

	if (node->child) {
		pll_rnetwork_node_t * snode = node->child;
		pll_rnetwork_node_t * new_snode = new_node;
		new_snode->child = (pll_rnetwork_node_t *) malloc(sizeof(pll_rnetwork_node_t));
		memcpy(new_snode->child, snode, sizeof(pll_rnetwork_node_t));
		new_snode->child->label = new_node->label;
		new_snode->child->reticulation_name = new_node->reticulation_name;
		snode = snode->child;
		new_snode = new_snode->child;
		new_snode->child = new_node;
	}
	if (node->left) {
		pll_rnetwork_node_t * snode = node->left;
		pll_rnetwork_node_t * new_snode = new_node;
		new_snode->left = (pll_rnetwork_node_t *) malloc(sizeof(pll_rnetwork_node_t));
		memcpy(new_snode->left, snode, sizeof(pll_rnetwork_node_t));
		new_snode->left->label = new_node->label;
		new_snode->left->reticulation_name = new_node->reticulation_name;
		snode = snode->left;
		new_snode = new_snode->left;
		new_snode->left = new_node;
	}
	if (node->right) {
		pll_rnetwork_node_t * snode = node->right;
		pll_rnetwork_node_t * new_snode = new_node;
		new_snode->right = (pll_rnetwork_node_t *) malloc(sizeof(pll_rnetwork_node_t));
		memcpy(new_snode->right, snode, sizeof(pll_rnetwork_node_t));
		new_snode->right->label = new_node->label;
		new_snode->right->reticulation_name = new_node->reticulation_name;
		snode = snode->right;
		new_snode = new_snode->right;
		new_snode->right = new_node;
	}

	return new_node;
}

static void rnetwork_recurse_clone(pll_rnetwork_node_t * new_root, const pll_rnetwork_node_t * root) {
	if (root->is_reticulation) {
		const pll_rnetwork_node_t * node = root->child;
		if (node) {
			new_root->child = clone_node(node);
			if (node->first_parent == root) {
				new_root->child->first_parent = new_root;
			} else {
				new_root->child->second_parent = new_root;
			}

			if (node->is_reticulation) {
				if (node->child) {
					pll_rnetwork_node_t * snode = node->child;
					pll_rnetwork_node_t * new_snode = new_root->child->child;
					rnetwork_recurse_clone(new_snode, snode);
				}
			} else {
				if (node->left) {
					pll_rnetwork_node_t * snode = node->left;
					pll_rnetwork_node_t * new_snode = new_root->child->left;
					rnetwork_recurse_clone(new_snode, snode);
				}
				if (node->right) {
					pll_rnetwork_node_t * snode = node->right;
					pll_rnetwork_node_t * new_snode = new_root->child->right;
					rnetwork_recurse_clone(new_snode, snode);
				}
			}
		}
	} else {
		const pll_rnetwork_node_t * node_left = root->left;
		if (node_left) {
			new_root->left = clone_node(node_left);
			new_root->left->parent = new_root;

			if (node_left->is_reticulation) {
				if (node_left->child) {
					pll_rnetwork_node_t * snode = node_left->child;
					pll_rnetwork_node_t * new_snode = new_root->left->child;
					rnetwork_recurse_clone(new_snode, snode);
				}
			} else {
				if (node_left->left) {
					pll_rnetwork_node_t * snode = node_left->left;
					pll_rnetwork_node_t * new_snode = new_root->left->left;
					rnetwork_recurse_clone(new_snode, snode);
				}
				if (node_left->right) {
					pll_rnetwork_node_t * snode = node_left->right;
					pll_rnetwork_node_t * new_snode = new_root->left->right;
					rnetwork_recurse_clone(new_snode, snode);
				}
			}
		}

		const pll_rnetwork_node_t * node_right = root->right;
		if (node_right) {
			new_root->right = clone_node(node_right);
			new_root->right->parent = new_root;

			if (node_right->is_reticulation) {
				if (node_right->child) {
					pll_rnetwork_node_t * snode = node_right->child;
					pll_rnetwork_node_t * new_snode = new_root->right->child;
					rnetwork_recurse_clone(new_snode, snode);
				}
			} else {
				if (node_right->left) {
					pll_rnetwork_node_t * snode = node_right->left;
					pll_rnetwork_node_t * new_snode = new_root->right->left;
					rnetwork_recurse_clone(new_snode, snode);
				}
				if (node_right->right) {
					pll_rnetwork_node_t * snode = node_right->right;
					pll_rnetwork_node_t * new_snode = new_root->right->right;
					rnetwork_recurse_clone(new_snode, snode);
				}
			}
		}
	}
}

PLL_EXPORT pll_rnetwork_node_t * pll_rnetwork_graph_clone(const pll_rnetwork_node_t * root) {
	pll_rnetwork_node_t * new_root = clone_node(root);
	const pll_rnetwork_node_t * snode = root;
	pll_rnetwork_node_t * new_snode = new_root;
	rnetwork_recurse_clone(new_snode, snode);
	return new_root;
}

PLL_EXPORT pll_rnetwork_t * pll_rnetwork_clone(const pll_rnetwork_t * network) {
	/* choose the last inner node as the starting point of the clone. It does not
	 really matter which node to choose, but since the newick parser places the
	 root node at the end of the list, we use the same notation here */
	pll_rnetwork_node_t * root = pll_rnetwork_graph_clone(network->root);
	return pll_rnetwork_wrapnetwork(root);
}
