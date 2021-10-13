program Main

  character(*), parameter :: InputFile = 'input.txt', OutputFile = 'data.plt' ! names of input and output files
  character MeshFile * 30        ! name of file with computational mesh
  integer, parameter :: IO = 12 ! input-output unit
  real, allocatable, dimension(:,:) :: X,Y,P,CellVolume ! scalar arrays
  real, allocatable, dimension(:,:,:) :: GradP, gradPExact, gradPError
  real, allocatable, dimension(:,:,:) :: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector ! vector arrays

!===  READ INPUT FILE ===
  write(*,*) 'Read input file: ', InputFile
  open(IO, FILE = InputFile)
  read(IO, *) MeshFile  ! read name of file with computational mesh
  close(IO)

!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
  write(*,*) 'Read nodes number from file: ', MeshFile
  open(IO, file = MeshFile)
  read(IO, *) NI, NJ
  write(*,*) 'NI, NJ = ', NI, NJ

!=== ALLOCATE ALL ARRAYS ===
  write(*,*) 'Allocate arrays'       
  allocate(X(NI,NJ)) ! mesh nodes X-coordinates
  allocate(Y(NI,NJ)) ! mesh nodes Y-coordinates
  allocate(P(0:NI,0:NJ))   ! Pressure
  allocate(GradP(0:NI, 0:NJ, 2)) ! Gradient of Pressure
  allocate(gradPExact(0:NI, 0:NJ, 2)) ! Exact value of gradient 
  allocate(gradPError(0:NI, 0:NJ, 2))
  allocate(CellVolume(NI-1,NJ-1))   ! Cell Volumes    
  allocate(CellCenter(0:NI,0:NJ,2)) ! Cell Centers
  allocate(IFaceCenter(NI,NJ-1,2)) ! Face Centers for I-faces
  allocate(IFaceVector(NI,NJ-1,2)) ! Face Vectors for I-faces
  allocate(JFaceCenter(NI-1,NJ,2)) ! Face Centers for J-faces
  allocate(JFaceVector(NI-1,NJ,2)) ! Face Vectors for J-faces

!===  READ GRID ===
  write(*,*) 'Read mesh from file: ', MeshFile
  read(IO,*) ((X(I,J), Y(I,J), I = 1, NI), J = 1, NJ)
  close(IO)

!=== CALCULATE METRIC ===
  write(*,*) 'Calculate metric'       
  call B_CalcMetric(NI, NJ, X, Y, CellCenter, CellVolume, IFaceCenter, IFaceVector, JFaceCenter, JFaceVector) 
  
!=== INITIATE FIELDS ===
  write(*,*) 'Initiate fields'       
  do  J = 0, NJ
    do  I = 0, NI
      P(I,J) = Pressure(CellCenter(I, J, 1), CellCenter(I, J, 2))
	  call CalcGradPExact(ni, nj, CellCenter(i, j, 1), CellCenter(i, j, 2), gradPExact)
    end do
  end do

!=== CALCULATE GRADIENT ===
  write(*,*) 'Calculate derivatives'    
  call B_CalcGradient(NI, NJ, P, GradP, CellVolume, CellCenter, IFaceCenter, JFaceCenter, IFaceVector, JFaceVector)
  gradPError = abs(GradP - gradPExact) / maxval(gradPExact)

!=== OUTPUT FIELDS ===
  write(*,*) 'Output fields to file: ', OutputFile       
  open(IO, file = OutputFile)
  call B_OutputFields(IO, NI, NJ, X, Y, P, GradP, gradPError)
  close(IO)

end program Main