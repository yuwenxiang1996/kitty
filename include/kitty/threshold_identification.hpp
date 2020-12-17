 /* kitty: C++ truth table library
 */

#pragma once

#include <vector>
#include <lpsolve/lp_lib.h> /* uncomment this line to include lp_solve */
#include "traits.hpp"

#include <iostream> // need to be deleted
using namespace std;

namespace kitty
{

/*! \brief Threshold logic function identification

  Given a truth table, this function determines whether it is a threshold logic function (TF)
  and finds a linear form if it is. A Boolean function is a TF if it can be expressed as

  f(x_1, ..., x_n) = \sum_{i=1}^n w_i x_i >= T

  where w_i are the weight values and T is the threshold value.
  The linear form of a TF is the vector [w_1, ..., w_n; T].

  \param tt The truth table
  \param plf Pointer to a vector that will hold a linear form of `tt` if it is a TF.
             The linear form has `tt.num_vars()` weight values and the threshold value
             in the end.
  \return `true` if `tt` is a TF; `false` if `tt` is a non-TF.
*/
template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
bool is_threshold( const TT& tt, std::vector<int64_t>* plf = nullptr )
{
  std::vector<int64_t> linear_form;
    
  // Get some paramters
  int32_t num_bits = static_cast<uint32_t>( tt.num_bits() );
  int32_t num_vars = static_cast<uint32_t>( tt.num_vars() );
  auto tt_flip = tt;
  int flip_table[50];
  int i, j, k, v, wk;
  // Print truth table
  std::cout << "Truth table of f is:";
  for ( int32_t i = num_bits - 1; i >= 0; i-- )
  {
    std::cout << get_bit( tt, i );
  }
  std::cout << endl;
  
// check unateness, create positive unate function
    for (auto i = 0u; i < num_vars; i++)
    {
        auto const tt_1 = cofactor0( tt_flip, i);
        auto const tt_2 = cofactor1( tt_flip, i);
        auto const s_t = tt_1 | tt_2;
        
        if (tt_1 = tt_2)
        {
            continue;
        }
        else if (tt_1 == s_t)
        {
            flip_inplace( tt_flip, i);
            neg_unate.push_back( i );
        }
        else if ( tt_2 != s_t)
        {
            return false;
        }
    }
  // Get cofactor
  //for ( i = 0u; i < num_vars ; ++i )
  //{
  //  auto fxi     = cofactor1(tt, i);
  //  auto not_fxi = cofactor0(tt, i);
  
    // Judge f is positive or negative unate
  //   if ( implies( fxi, not_fxi ) )
  //  {
  //    std::cout << "f is negative unate in x" << i << endl;
  //    tt_flip = flip( tt_flip, i ); // Flip xi so f_flip is positive unate in xi
  //   flip_table[i] = 1;
  //  }
  //  else if ( implies( not_fxi, fxi ) )
  //  {
  //    std::cout << "f is positive unate in x" << i << endl;
  //    flip_table[i] = 0;
   // }
  //  else
  //  {
  //    std::cout << "f is binate in x" << i << endl;
  //    return false;
  //  }
 // }
  // check tt_flip
    
  // Print truth table of f_flip
  //std::cout << "Truth table of f_flip is:";
  //for (  i = num_bits - 1; i >= 0; i-- )
  //{
  //  std::cout << get_bit( tt_flip, i );
  //}
  //std::cout << endl;

  /* Create a New ILP problem with tt_flip*/
  lprec* lp;
  int Ncol, *colno = NULL, ret = 0;
  REAL* row = NULL;

  Ncol = num_vars + 1;
  lp = make_lp( 0, Ncol );

  // Set var names
  if ( ret == 0 )
  {
    for ( k = 0; k < num_vars; k++ )
    {
      set_col_name( lp, k+1, "C1" );
    }
    set_col_name( lp, num_vars + 1, "T" );

    // create space large enough for one row 
    colno = (int*)malloc( Ncol * sizeof( *colno ) );
    row = (REAL*)malloc( Ncol * sizeof( *row ) );
    if ( ( colno == NULL ) || ( row == NULL ) )
      ret = 2;
  }
  for ( j = 0; j < num_bits; j++ )
  {
    v = j;
    for ( k = 0; k < num_vars; k++ )
    {
      colno[k] = k + 1; /* k column */
      
      wk = v % 2;
      v = v >> 1;
      //std::cout << wk;
      row[k] = wk? 1 : 0;
    }
    //std::cout << endl;
    col[num_vars] = num_vars + 1; /* num_vars column */
    row[num_vars] = -1;
    if ( get_bit(tt_flip,j) )
        add_constraintex( lp, num_vars + 1, row, colno, GE, 0 );
    else
        add_constraintex( lp, num_vars + 1, row, colno, LE, -1 );
  }
  // Set the objective function
  if ( ret == 0 )
  {
    set_add_rowmode( lp, FALSE ); /* rowmode should be turned off again when finished building the model */

    for ( k = 0; k <= num_vars; k++ )
    {
      colno[k] = k + 1; /* k column */
      row[k] = 1;
    }
    /* set the objective in lpsolve */
    if ( !set_obj_fnex( lp, num_vars + 1, row, colno ) )
      ret = 4;
  }

  // Set variables to int
  if ( ret == 0 )
  {
    for ( k = 1; k <= num_vars + 1; k++ )
      set_int( lp, k, TRUE );
  }

  if ( ret == 0 )
  {
    /* set the object direction to minimize */
    set_minim( lp );

    /* just out of curioucity, now show the model in lp format on screen */
    write_LP( lp, stdout );

    /* I only want to see important messages on screen while solving */
    set_verbose( lp, IMPORTANT );

    /* Now let lpsolve calculate a solution */
    ret = solve( lp );
    if ( ret == OPTIMAL )
      ret = 0;
    else
      ret = 5;
  }

  if ( ret == 0 )
  {
    /* objective value */
    printf( "Objective value: %f\n", get_objective( lp ) );

    /* variable values */
    get_variables( lp, row );
    int T = row[Ncol - 1];
    for ( j = 0; j < Ncol; j++ )
    {
      printf( "%s: %f\n", get_col_name( lp, j + 1 ), row[j] );
    }
    for ( j = 0; j < num_vars; j++ )
    {
      if ( flip_table[j] == 0 )
      {
        linear_form.push_back( row[j] );
      }
      else
      {
        linear_form.push_back( -row[j] );
        T = T - row[j];
      }
    }
    linear_form.push_back( T );
  }

  if ( row != NULL )
    free( row );
  if ( col != NULL )
    free( colno );

  /* clean up such that all used memory by lpsolve is free */
  if ( lp != NULL )
    delete_lp( lp );

  /* if tt is non-TF: */
  if (ret != 0 )
    return false;

  /* if tt is TF: */
  /* push the weight and threshold values into `linear_form` */
  if ( plf )
  {
    *plf = linear_form;
  }
  return true;
}

} /* namespace kitty */
