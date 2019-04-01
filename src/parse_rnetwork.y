/*
    Copyright (C) 2015 Tomas Flouri, with changes and additions by Sarah Lutteropp in 2019

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
%{
#include "pll.h"

extern int pll_rnetwork_lex();
extern FILE * pll_rnetwork_in;
extern void pll_rnetwork_lex_destroy();
extern int pll_rnetwork_lineno;
extern int pll_rnetwork_colstart;
extern int pll_rnetwork_colend;

extern int pll_rnetwork_parse();
extern struct pll_rnetwork_buffer_state * pll_rnetwork__scan_string(const char * str);
extern void pll_rnetwork__delete_buffer(struct pll_rnetwork_buffer_state * buffer);

static unsigned int tip_cnt = 0;
static unsigned int inner_tree_cnt = 0;
static unsigned int reticulation_cnt = 0;
static pll_rnetwork_node_t ** reticulation_node_pointers;
static char* reticulation_node_names[64];

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

    free(node);
  }

  /* deallocate network structure */
  free(network->nodes);
  free(network);
}

static void pll_rnetwork_error(pll_rnetwork_node_t * node, const char * s)
{
  pll_errno = PLL_ERROR_NEWICK_SYNTAX;
  if (pll_rnetwork_colstart == pll_rnetwork_colend)
    snprintf(pll_errmsg, 200, "%s. (line %d column %d)\n",
             s, pll_rnetwork_lineno, pll_rnetwork_colstart);
  else
    snprintf(pll_errmsg, 200, "%s. (line %d column %d-%d)\n",
             s, pll_rnetwork_lineno, pll_rnetwork_colstart, pll_rnetwork_colend);
}

%}

%union
{
  char * s;
  char * d;
  struct pll_rnetwork_node_s * network;
}

%error-verbose
%parse-param {struct pll_rnetwork_node_s * network}
%destructor { pll_rnetwork_graph_destroy($$,NULL); } subtree
%destructor { free($$); } STRING
%destructor { free($$); } NUMBER
%destructor { free($$); } label

%token<s> STRING
%token<d> NUMBER
%type<s> label optional_label
%type<d> number optional_length optional_number_after_colon
%type<network> subtree
%start input
%%

input: '(' subtree ',' subtree ')' optional_label optional_length ';'
{
  inner_tree_cnt++;
  network->is_reticulation = 0;

  network->left   = $2;
  network->right  = $4;
  network->label  = $6;
  network->length = $7 ? atof($7) : 0;
  network->parent = NULL;
  free($7);

  if (network->left->is_reticulation)
  {
    if (!network->left->first_parent)
    {
      network->left->first_parent = network;
    }
    else
    {
      network->left->second_parent = network;
    }
  }
  else
  {
    network->left->parent = network;
  }
  
  if (network->right->is_reticulation)
  {
    if (!network->right->first_parent)
    {
      network->right->first_parent = network;
    }
    else
    {
      network->right->second_parent = network;
    }
  }
  else
  {
    network->right->parent = network;
  }
};

subtree: '(' subtree ',' subtree ')' optional_label optional_length
{
  inner_tree_cnt++;
  $$ = (pll_rnetwork_node_t *)calloc(1, sizeof(pll_rnetwork_node_t));
  $$->is_reticulation = 0;
  $$->left   = $2;
  $$->right  = $4;
  $$->label  = $6;
  $$->length = $7 ? atof($7) : 0;
  free($7);

  $$->left->parent  = $$;
  $$->right->parent = $$;

  if ($$->left->is_reticulation)
  {
    if (!$$->left->first_parent)
    {
      $$->left->first_parent = $$;
    }
    else
    {
      $$->left->second_parent = $$;
    }
  }
  else
  {
    $$->left->parent = $$;
  }
  
  if ($$->right->is_reticulation)
  {
    if (!$$->right->first_parent)
    {
      $$->right->first_parent = $$;
    }
    else
    {
      $$->right->second_parent = $$;
    }
  }
  else
  {
    $$->right->parent = $$;
  }


}
       | '(' subtree ')' optional_label '#' label ':' optional_number_after_colon ':' optional_number_after_colon ':' number
{
  $$ = (pll_rnetwork_node_t *)calloc(1, sizeof(pll_rnetwork_node_t));
  $$->is_reticulation = 1;
  $$->child   = $2;
  $$->left   = NULL;
  $$->right  = NULL;
  $$->first_parent   = NULL;
  $$->second_parent  = NULL;
  $$->label  = $4;
  $$->reticulation_name = $6;
  $$->reticulation_index = reticulation_cnt;
  $$->length = $8 ? atof($8) : 0;
  free($8);
  $$->support = $10 ? atof($10) : 0;
  free($10);
  $$->prob = atof($12);
  free($12);

  $$->child->parent  = $$;

  reticulation_node_pointers[reticulation_cnt] = $$;
  reticulation_node_names[reticulation_cnt] = $$->reticulation_name;
  reticulation_cnt++;
}
       | '(' subtree ')' optional_label '#' label
{
  $$ = (pll_rnetwork_node_t *)calloc(1, sizeof(pll_rnetwork_node_t));
  $$->is_reticulation = 1;
  $$->child   = $2;
  $$->left   = NULL;
  $$->right  = NULL;
  $$->first_parent   = NULL;
  $$->second_parent  = NULL;
  $$->label  = $4;
  $$->reticulation_name = $6;
  $$->reticulation_index = reticulation_cnt;
  $$->length = 0;
  $$->support = 0;
  $$->prob = 0;

  $$->child->parent  = $$;
  reticulation_node_pointers[reticulation_cnt] = $$;
  reticulation_node_names[reticulation_cnt] = $$->reticulation_name;
  reticulation_cnt++;
}
       | optional_label '#' label optional_length
{
  unsigned int i = 0;
  for (i = 0; i < reticulation_cnt; ++i)
  {
    if (strcmp(reticulation_node_names[i],$3) == 0)
    {
      $$ = reticulation_node_pointers[i];
      break;
    }
  }
  free($4);
}
       | label optional_length
{
  $$ = (pll_rnetwork_node_t *)calloc(1, sizeof(pll_rnetwork_node_t));
  $$->is_reticulation = 0;
  $$->label  = $1;
  $$->length = $2 ? atof($2) : 0;
  $$->left   = NULL;
  $$->right  = NULL;
  $$->first_parent   = NULL;
  $$->second_parent  = NULL;
  $$->child = NULL;
  tip_cnt++;
  free($2);
};


optional_label:  {$$ = NULL;} | label  {$$ = $1;};
optional_length: {$$ = NULL;} | ':' number {$$ = $2;};
label: STRING    {$$=$1;} | NUMBER {$$=$1;};
number: NUMBER   {$$=$1;};

optional_number_after_colon: {$$ = NULL;} | number {$$ = $1;};

%%

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
                                 unsigned int * tip_index,
                                 unsigned int * inner_index)
{
  if (!node)
  {
    return;
  }
  
  if (!node->is_reticulation && !node->left && !node->right)
  { // we are at a tip node
    node->idx = *tip_index;
    array[*tip_index] = node;
    *tip_index = *tip_index + 1;
    return;
  }
  
  if (!node->is_reticulation)
  {
    fill_nodes_recursive(node->left,  array, tip_index, inner_index);
    fill_nodes_recursive(node->right, array, tip_index, inner_index);
  }
  else
  {
    fill_nodes_recursive(node->child, array, tip_index, inner_index);
  }

  array[*inner_index] = node;
  node->idx = *inner_index;
  *inner_index = *inner_index + 1;
}

/*static unsigned int rnetwork_count_tips(pll_rnetwork_node_t * root)
{
  unsigned int count = 0;

  if (root->left)
    count += rtree_count_tips(root->left);
  if (root->right)
    count += rtree_count_tips(root->right);

  if (!root->left && !root->right)
    return 1;

  return count;
}*/

PLL_EXPORT pll_rnetwork_t * pll_rnetwork_wrapnetwork(pll_rnetwork_node_t * root)
{
  pll_rnetwork_t * network = (pll_rnetwork_t *)malloc(sizeof(pll_rnetwork_t));
  if (!network)
  {
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    pll_errno = PLL_ERROR_MEM_ALLOC;
    return PLL_FAILURE;
  }
  
  network->tip_count = tip_cnt;
  network->inner_tree_count = inner_tree_cnt;
  network->reticulation_count = reticulation_cnt;

  //if (tip_count < 2 && tip_count != 0)
  //{
    //snprintf(pll_errmsg, 200, "Invalid tip_count value (%u).", tip_count);
    //pll_errno = PLL_ERROR_PARAM_INVALID;
    //return PLL_FAILURE;
  //}

  //if (tip_count == 0)
  //{
    //tip_count = rnetwork_count_tips(root);
    //if (tip_count < 2)
    //{
      //snprintf(pll_errmsg, 200, "Input tree contains no inner nodes.");
      //pll_errno = PLL_ERROR_PARAM_INVALID;
      //return PLL_FAILURE;
    //}
  //}

  unsigned int total_node_cnt = tip_cnt + inner_tree_cnt + reticulation_cnt;
  network->nodes = (pll_rnetwork_node_t **)malloc(total_node_cnt * sizeof(pll_rnetwork_node_t *));
  if (!network->nodes)
  {
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    pll_errno = PLL_ERROR_MEM_ALLOC;
    return PLL_FAILURE;
  }
  
  unsigned int tip_index = 0;
  unsigned int inner_index = tip_cnt;

  fill_nodes_recursive(root->left, network->nodes, &tip_index, &inner_index);
  fill_nodes_recursive(root->right, network->nodes, &tip_index, &inner_index);
  network->nodes[inner_index] = root;
  network->root = root;

  return network;
}

PLL_EXPORT pll_rnetwork_t * pll_rnetwork_parse_newick(const char * filename)
{
  pll_rnetwork_t * network;

  struct pll_rnetwork_node_s * root;

  /* reset counters */
  tip_cnt = 0;
  inner_tree_cnt = 0;
  reticulation_cnt = 0;
  reticulation_node_pointers = (pll_rnetwork_node_t **)calloc(64, sizeof(pll_rnetwork_node_t*));

  /* open input file */
  pll_rnetwork_in = fopen(filename, "r");
  if (!pll_rnetwork_in)
  {
    pll_errno = PLL_ERROR_FILE_OPEN;
    snprintf(pll_errmsg, 200, "Unable to open file (%s)", filename);
    
    /* free the counters */
    free(reticulation_node_pointers);
    
    return PLL_FAILURE;
  }

  /* create root node */
  if (!(root = (pll_rnetwork_node_t *)calloc(1, sizeof(pll_rnetwork_node_t))))
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    
    /* free the counters */
    free(reticulation_node_pointers);
    
    return PLL_FAILURE;
  }

  if (pll_rnetwork_parse(root))
  {
    pll_rnetwork_graph_destroy(root,NULL);
    root = NULL;
    fclose(pll_rnetwork_in);
    pll_rnetwork_lex_destroy();
    
    /* free the counters */
    free(reticulation_node_pointers);
    
    return PLL_FAILURE;
  }

  if (pll_rnetwork_in) fclose(pll_rnetwork_in);

  pll_rnetwork_lex_destroy();

  //initialize clv and scaler indices
  //pll_rnetwork_reset_template_indices(root, tip_cnt);

  /* wrap network */
  network = pll_rnetwork_wrapnetwork(root);
  
  /* free the counters */
  free(reticulation_node_pointers);

  return network;
}

PLL_EXPORT pll_rnetwork_t * pll_rnetwork_parse_newick_string(const char * s)
{
  int rc;
  struct pll_rnetwork_node_s * root;
  pll_rnetwork_t * network = NULL;

  /* reset counters */
  tip_cnt = 0;
  inner_tree_cnt = 0;
  reticulation_cnt = 0;
  reticulation_node_pointers = (pll_rnetwork_node_t **)calloc(64, sizeof(pll_rnetwork_node_t*));

  if (!(root = (pll_rnetwork_node_t *)calloc(1, sizeof(pll_rnetwork_node_t))))
  {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "Unable to allocate enough memory.");
    
    /* free the counters */
    free(reticulation_node_pointers);
    
    return PLL_FAILURE;
  }

  struct pll_rnetwork_buffer_state * buffer = pll_rnetwork__scan_string(s);
  rc = pll_rnetwork_parse(root);
  pll_rnetwork__delete_buffer(buffer);

  pll_rnetwork_lex_destroy();

  if (!rc)
  {
    //initialize clv and scaler indices */
    //pll_rnetwork_reset_template_indices(root, tip_cnt);

    network = pll_rnetwork_wrapnetwork(root);
  }
  else
    free(root);
    
  /* free the counters */
  free(reticulation_node_pointers);

  return network;
}
