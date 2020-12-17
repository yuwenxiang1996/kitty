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
template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
bool is_threshold( const TT& tt, std::vector<int64_t>* plf = nullptr )
{
  std::vector<int64_t> linear_form;
    
  // Get some paramters
  int32_t num_bits = static_cast<uint32_t>( tt.num_bits() );
  int32_t num_vars = static_cast<uint32_t>( tt.num_vars() );
  
  auto ttf = tt;
  int flip_table[50];
  int i, j, m, v, wk;
    
    
  // Print truth table
  /*std::cout << "Truth table of f is:";
  for ( int32_t i = num_bits - 1; i >= 0; i-- )
  {
    std::cout << get_bit( tt, i );
  }
  std::cout << endl;
 */
 
    //check unateness, create positive unate function
    for (auto i = 0u; i < num_vars; i++)
    {
        auto const tt_1 = cofactor0( ttf, i);
        auto const tt_2 = cofactor1( ttf, i);
        auto const s_t = tt_1 | tt_2;
        
        if ( implies( tt_2 , tt_1 ) )  // if positive unate
        {
            ttf = flip ( ttf, i);
            flip_table[i] = 1;
        }
        else if ( implies (tt_1 , tt_2 ) ) // if negative unate
        {
            flip_table[i] = 0;
        }
        else if (tt_2 != s_t)         // if binate
        {
            /* binate tt is non-TF*/
            return false;
        }
    }
 
  // Get cofactor
  /*for ( i = 0u; i < num_vars ; ++i )
  {
    auto fxi     = cofactor1(tt, i);
    auto not_fxi = cofactor0(tt, i);
  
    // Judge f is positive or negative unate
     if ( implies( fxi, not_fxi ) )
    {
      std::cout << "f is negative unate in x" << i << endl;
      tt_flip = flip( tt_flip, i ); // Flip xi so f_flip is positive unate in xi
     flip_table[i] = 1;
    }
    else if ( implies( not_fxi, fxi ) )
    {
      std::cout << "f is positive unate in x" << i << endl;
      flip_table[i] = 0;
    }
    else
    {
      std::cout << "f is binate in x" << i << endl;
      return false;
    }
  }
  // check tt_flip
  */
   //Print truth table of f_flip
  /*std::cout << "Truth table of f_flip is:";
  for (  i = num_bits - 1; i >= 0; i-- )
  {
    std::cout << get_bit( tt_flip, i );
  }
  std::cout << endl;
   

  /* Create a New ILP problem with tt_flip */
  lprec* lp;
  int Ncol, *colno = NULL, ret = 0;
  REAL* row = NULL;

  Ncol = num_vars + 1;
  lp = make_lp( 0, Ncol );

  // Set var names
  if ( ret == 0 )
  {
    for ( m = 0; m < num_vars; m++ )
    {
      set_col_name( lp, m+1, "C1" );
    }
    set_col_name( lp, num_vars + 1, "T" );

    // create space large enough for one row 
    colno = (int*)malloc( Ncol * sizeof( *colno ) );
    row = (REAL*)malloc( Ncol * sizeof( *row ) );
    if ( ( colno == NULL ) || ( row == NULL ) )
      ret = 1;
  }
  for ( j = 0; j < num_bits; j++ )
  {
    v = j;
    for ( m = 0; m < num_vars; m++ )
    {
      colno[m] = m + 1; /* m column */
      
      wk = v % 2;
      v = v >> 1;
      //std::cout << wk;
      row[m] = wk? 1 : 0;
    }
    //std::cout << endl;
    colno[num_vars] = num_vars + 1; /* num_vars column */
    row[num_vars] = -1;
    if ( get_bit(ttf,j) )
        add_constraintex( lp, num_vars + 1, row, colno, GE, 0 );
    else
        add_constraintex( lp, num_vars + 1, row, colno, LE, -1 );
  }
  // Set the objective function
  if ( ret == 0 )
  {
    set_add_rowmode( lp, FALSE ); /* rowmode should be turned off again when finished building the model */

    for ( m = 0; m <= num_vars; m++ )
    {
      colno[m] = m + 1; /* k column */
      row[m] = 1;
    }
    /* set the objective in lpsolve */
    if ( !set_obj_fnex( lp, num_vars + 1, row, colno ) )
      ret = 2;
  }

  // Set variables to int
  if ( ret == 0 )
  {
    for ( m = 1; m <= num_vars + 1; m++ )
      set_int( lp, m, TRUE );
  }

  if ( ret == 0 )
  {
    set_minim( lp );

    write_LP( lp, stdout );

    set_verbose( lp, IMPORTANT );

    /*  lpsolve calculate a solution */
    ret = solve( lp );
    if ( ret == OPTIMAL )
      ret = 0;
    else
      ret = 3;
  }

  if ( ret == 0 )
  {
    /* objective value */
    //printf( "Objective value: %f\n", get_objective( lp ) );

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

 /* if tt is non-TF */
  if (ret != 0 )
  {
    if ( row != NULL )
     free( row );
    if ( colno != NULL )
     free( colno );
    if ( lp != NULL )
        delete_lp( lp);
    return false;
  }
    

  /* if tt is TF: */
  /* push the weight and threshold values into `linear_form` */
  if ( plf )
  {
    *plf = linear_form;
  }
  return true;
}

} /* namespace kitty */
