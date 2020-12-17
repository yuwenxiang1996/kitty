 /* kitty: C++ truth table library
 */

#pragma once

#include <vector>
#include <string>
#include <bitset>
#include <algorithm>
#include <lpsolve/lp_lib.h> /* uncomment this line to include lp_solve */
#include "traits.hpp"
#include "properties.hpp"
#include "dynamic_truth_table.hpp"
#include "static_truth_table.hpp"


//#include <iostream> // need to be deleted
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
  int i, j, k;
    
 
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
        else if ( tt_2 != s_t )      // if binate
        {
            /* binate tt is non-TF*/
            return false;
        }
    }

    /* build ILP problem */
    lprec* lp;
       int *colno = NULL, ret = 0;
       uint64_t Nrow;
       REAL* row = NULL;
    
    const uint64_t NCol = ttf.num_vars() + 1;
    Nrow = pow( 2, ttf.num_vars() );
    lp = make_lp(0, NCol );
    
    if (lp == NULL )
        ret = 1; // fail to construct
    
    if ( ret == 0 )
    {
        /* make sure a large enough space for row*/
        colno = (int*)malloc( NCol * sizeof( *colno ));
        row = (REAL*)malloc( NCol * sizeof( *row ));
        /* fail to construct */
        if ( ( colno == NULL ) || ( row == NULL ) )
        ret = 2;
    }
    
    if ( ret == 0 )
    {
        set_add_rowmode( lp, TRUE);
        
        for ( k = 0; k < Nrow; k++ )
        {
            j = 0 ;
            std::bitset<64> bs( k );
            for ( uint64_t var = 1; var <= NCol -1; var++ )
            {
                colno[j] = var;
                if ( bs.test( var - 1 ))
                {
                    row[j++] = 1;
                }
                else
                {
                    row[j++] = 0;
                }
            }
            
            colno[j] = NCol; /* last column*/
            row[j++] = -1;
            
            if ( get_bit( ttf, k ))
            {
                if ( !add_constraintex( lp, j, row, colno, GE, 0))
                    ret = 3;
            }
            else
            {
                if ( !add_constraintex( lp, j, row, colno, LE, -1))
                    ret = 3;
            }
        }
    }
    
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
      for ( k = 1; k <= num_vars + 1; k++ )
    {
          set_int( lp, k, TRUE );
    }

    if ( ret == 0 )
    {
      set_minim( lp );

      set_verbose( lp, IMPORTANT );

      /*  lpsolve calculate a solution */
      ret = solve( lp );
        
      //res = ret;
        
      if ( ret == OPTIMAL )
        ret = 0;
      else
        ret = 5;
    }
    
    if ( ret == 0 )
    {
        get_variables( lp, row );
        int T = row[NCol - 1];
        for ( j = 0 ; j < num_vars; j++)
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
    /* free allocated memory*/
    if ( row != NULL )
     free( row );
    if ( colno != NULL )
     free( colno );
    /* free all used memory used by lpsolve */
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

}
/* namespace kitty */
