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
%destructor { pll_rnetwork_graph_destroy($$,NULL); } subnetwork
%destructor { free($$); } STRING
%destructor { free($$); } NUMBER
%destructor { free($$); } label

%token<s> STRING
%token<d> NUMBER
%type<s> label optional_label
%type<d> number optional_length optional_number_after_colon
%type<network> subnetwork
%start input
%%

input: '(' subnetwork ',' subnetwork ')' optional_label optional_length ';'
{
  inner_tree_cnt++;
  network->is_reticulation = 0;

  network->left   = $2;
  network->right  = $4;
  network->label  = $6;
  network->length = $7 ? atof($7) : 0;
  network->parent = NULL;
  network->idx = 0;
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

subnetwork: '(' subnetwork ',' subnetwork ')' optional_label optional_length
{
  inner_tree_cnt++;
  $$ = (pll_rnetwork_node_t *)calloc(1, sizeof(pll_rnetwork_node_t));
  $$->is_reticulation = 0;
  $$->left   = $2;
  $$->right  = $4;
  $$->label  = $6;
  $$->length = $7 ? atof($7) : 0;
  free($7);
  $$->idx = 0;

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
       | '(' subnetwork ')' optional_label '#' label ':' optional_number_after_colon ':' optional_number_after_colon ':' number
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
  $$->first_parent_length = $8 ? atof($8) : 0;
  free($8);
  $$->support = $10 ? atof($10) : 0;
  free($10);
  $$->prob = atof($12);
  free($12);
  $$->idx = 0;

  $$->child->parent  = $$;

  reticulation_node_pointers[reticulation_cnt] = $$;
  reticulation_node_names[reticulation_cnt] = $$->reticulation_name;
  reticulation_cnt++;
}
       | '(' subnetwork ')' optional_label '#' label
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
  $$->idx = 0;

  $$->child->parent  = $$;
  reticulation_node_pointers[reticulation_cnt] = $$;
  reticulation_node_names[reticulation_cnt] = $$->reticulation_name;
  reticulation_cnt++;
}
       | optional_label '#' label ':' optional_number_after_colon ':' optional_number_after_colon ':' number // branch length, support, probability
{
  unsigned int i = 0;
  for (i = 0; i < reticulation_cnt; ++i)
  {
    if (strcmp(reticulation_node_names[i],$3) == 0)
    {
      $$ = reticulation_node_pointers[i];
      if ($5) {
        $$->second_parent_length = atof($5);
        free($5);
      } else {
        $$->second_parent_length = 0;
      }
      if ($7) {
        $$->support = atof($7);
        free($7); 
      } else {
        $$->support = 0;
      }
      if ($9) {
        $$->prob = atof($9);
        free($9); 
      } else {
        $$->prob = 0.5;
      }
      break;
    }
  }
}
       | optional_label '#' label optional_length
{
  unsigned int i = 0;
  for (i = 0; i < reticulation_cnt; ++i)
  {
    if (strcmp(reticulation_node_names[i],$3) == 0)
    {
      $$ = reticulation_node_pointers[i];
      if ($4) {
        $$->second_parent_length = atof($4);
        free($4);
      } else {
        $$->second_parent_length = 0;
      }
      $$->support = 0;
      break;
    }
  }
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
  $$->idx = 0;
  tip_cnt++;
  free($2);
};


optional_label:  {$$ = NULL;} | label  {$$ = $1;};
optional_length: {$$ = NULL;} | ':' number {$$ = $2;};
label: STRING    {$$=$1;} | NUMBER {$$=$1;};
number: NUMBER   {$$=$1;};

optional_number_after_colon: {$$ = NULL;} | number {$$ = $1;};

%%

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
