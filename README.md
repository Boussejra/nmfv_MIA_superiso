# superiso_nmfv
Working version of code for nmfv version of superiso

## Dependencies 

This particular version of Superiso uses the gsl library for some computations. Make sure you have it isntalled and on Path. 

## Usage

To use the code immediately, you can launch `./cw_nmfv_slha.x <inputfile.lha>` with `<inputfile.lha>` your SLHA input file. You may need to recompile some files first.

The SLHA input file must further have an additional input block for the MIA parameters : 

```
BLOCK NMFVPAR   #NMFV params
 1      <val>     # delta_u_23_LR
 2      <val>     # delta_u_23_LL 
 3      <val>     # delta_u_33_LR
 4      <val>     # delta_d_23_LL 
 5      <val>     # delta_d_23_RR
 6      <val>     # delta_d_23_RL
 7      <val>     # delta_d_23_LR 
 8      <val>     # delta_d_33_RL 
 9      <val>     # delta_d_33_LR
```
With `<val>` a floating point value in the range [-1,1]
An example is provided in `example_MIA.lha`. 

## Output

The code outputs C7, C9, C10 both at the mu= m_W scale and at the mu=m_b scale.
The output stream is currently set to `out.txt`. 


