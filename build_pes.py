"""Extract unchanged scalar LSTH routines, not a refit. Upstream MIT header retained."""
from pathlib import Path
import hashlib, re, json
p=Path(__file__).parent
s=(p/'vendor/H3_m.f90').read_text()
header=s[:s.index('MODULE QML_H3_m')]
constants=s[s.index('    real(kind=Rkind), parameter :: C6'):s.index('!> @brief Derived type')]
routines=s[s.index('  SUBROUTINE QML_LSTH_refactoring'):s.index('  SUBROUTINE EvalFunc_QML_H3(')]
wrapper='''
subroutine potential(n,xyz,v,legacy) bind(C)
use iso_c_binding
integer(c_int),value::n,legacy
real(c_double),intent(in)::xyz(3,n)
real(c_double),intent(out)::v(n)
integer::i
do i=1,n
 if(legacy==0)then
 call QML_LSTH_refactoring(xyz(:,i),v(i))
 else
 call QML_LSTH(xyz(:,i),v(i))
 endif
enddo
end subroutine
'''
text=header+'module scalar_lsth\nimplicit none\ninteger,parameter::Rkind=kind(1d0)\nreal(Rkind),parameter::ZERO=0d0,ONE=1d0,TWO=2d0,SIX=6d0,EIGHT=8d0,TEN=10d0\n'+constants+'\ncontains\n'+routines+wrapper+'end module\n'
(p/'vendor/scalar_lsth.f90').write_text(text)
(p/'vendor/EXTRACTION.json').write_text(json.dumps({'source_commit':'73f2ef3d322b253a4ef39e08216261bd974180f6','upstream_sha256':hashlib.sha256(s.encode()).hexdigest(),'operation':'unchanged constants and scalar routines; iso_c_binding vector wrapper; no parameter changes'},indent=2))
