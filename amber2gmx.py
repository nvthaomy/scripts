import parmed as pmd

amber = pmd.load_file('500opc_1nacl.parm7','500opc_1nacl.crd')
amber.save('500opc_1nacl.top')
amber.save('500opc_1nacl.gro')

