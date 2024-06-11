module type_kinds
    use iso_fortran_env, only:  int8, int16, int32, int64, real32, real64, &
                                input_unit, output_unit, error_unit
    implicit none
    integer, parameter :: sp        = real32
    integer, parameter :: dp        = real64
    integer, parameter :: stdin     = input_unit
    integer, parameter :: stdout    = output_unit
    integer, parameter :: stderr    = error_unit
end module type_kinds