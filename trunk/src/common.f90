module commonFunctions

	type dstring
		character (len=80) :: data
		integer :: length
	end type dstring

	type pstring
		character, dimension(:), pointer :: data
		integer :: length
	end type pstring

	interface new
		module procedure dstring_init	
			end interface

    interface newp
        module procedure pstring_init
    end interface

    !interface reallocate
    !   module procedure reallocate_pstr    	
    !end interface

    interface suggestFileName
        module procedure suggestFileName_str, suggestFileName_dstr
    end interface

    interface lmod
        module procedure lmod_int, lmod_real4, lmod_real8
    end interface

    interface replace
        module procedure replace_str, replace_dstr
    end interface

    interface operator(+)
        module procedure dstring_concat
    end interface

    contains

         function pstring_init(str)
            implicit none

            ! Function type
             type(pstring) pstring_init

            ! Args
             character(len=*) str

            ! Locals
             integer istrLen1, istrLen2

             iStrLen1 = len(str)
             iStrLen2 = strLen(str)

														if(iStrLen1 > 0 .and. iStrLen2 == 0) then
             ! Probably just spaces
															allocate(pstring_init%data(iStrLen1))
															pstring_init%length = iStrLen1
															pstring_init%data(1:iStrLen1) = str(1:iStrLen1)

														else
             ! Enter data from str
															allocate(pstring_init%data(iStrLen2))
															pstring_init%length = iStrLen2
															pstring_init%data(1:iStrLen2) = str(1:iStrLen2)

														endif

														return

													end function pstring_init	


													function dstring_init(str)
														implicit none

          ! Function type
														type(dstring) dstring_init

          ! Args
														character(len=*) str

          ! Locals
														integer istrLen1, istrLen2

														iStrLen1 = len(str)
														iStrLen2 = strLen(str)

														if(iStrLen1 > 0 .and. iStrLen2 == 0) then
             ! Probably just spaces
															dstring_init%length = iStrLen1
															dstring_init%data = str(1:iStrLen1)

														else
             ! Enter data from str
															dstring_init%length = iStrLen2
															dstring_init%data = str(1:iStrLen2)

														endif

		! eliminate any remaining stuff on the memory
														if(dstring_init%length < 80) then
															dstring_init%data(dstring_init%length+1:80) = ''
														endif
														return
													end function dstring_init

													function dstring_concat(dstr1,dstr2)																
														implicit none

         ! Function delcaration
														type(dstring) dstring_concat

         ! Args
														type(dstring) , intent(in) :: dstr1,dstr2

														dstring_concat%data = dstr1%data(1:dstr1%length) // dstr2%data(1:dstr2%length)
														dstring_concat%length = dstr1%length + dstr2%length

														return																
													end function dstring_concat



													function dStrToReal(dStr)
														implicit none

          ! Function type
														real dStrToReal

          ! Args
														type(dstring) dStr

          ! Locals
														integer iTemp
														real rTemp

														rTemp = 9999999.9999
														iTemp = getUnit()
														open(unit=iTemp,status='SCRATCH')
														write(iTemp,*) dStr%data(1:dStr%length)
														rewind(iTemp)
														read(iTemp, *) rTemp
														close(iTemp)

														if(abs(rTemp - 9999999.9999) > 0.000001) then
															dStrToReal = rTemp
														else
															call warn("dStrToReal()", " Could not assign a value to the real number.")		 	 
														end if

														return

													end function dStrToReal


													function dStrToStr(dStr)
														implicit none

          ! Function type
														character (len=80) dStrToStr

          ! Args
														type(dstring) dStr

          ! Locals
														integer i

														dStrToStr(1:80) = ''

														do i = 1, dStr%length
															dStrToStr(i:i) = dStr%data(i:i)
														enddo

														return

													end function dStrToStr

													function dStrToInt(dStr)
														implicit none

          ! Function type
														integer dStrToInt

          ! Args
														type(dstring) dStr

          ! Locals
														integer iTemp
														integer iInt


														iTemp = getUnit()
														open(unit=iTemp,status='SCRATCH')
														write(iTemp,*) dStr%data(1:dStr%length)
														rewind(iTemp)
														read(iTemp, *) iInt
														close(iTemp)


														dStrToInt = iInt

														return

													end function dStrToInt

													subroutine dStrWrite(iunit,dstr)
														implicit none

			!Args
														integer :: iunit
														type(dstring) :: dstr

														write(iunit,*) dstr%data(1:dstr%length)

														return
													end subroutine dStrWrite

													subroutine createStrArr(str,dArr,iNum,iCount)
														implicit none

			! Args
														character (len=*) str
														integer iNum
														type(dstring) dArr(iNum)
														integer, optional, intent(out) :: iCount

			! Locals
														integer i,j
														integer iChar
														logical foundChar
														type(dstring) dStr, dTmp


			! Set up dstring
														dStr = new(str)

			! Set initial statements
														foundChar = .false.
														iChar = 0
														j = 1

			! Loop over the string and grab anything that is there
														do i = 1, dStr%length
															if(dStr%data(i:i) .ne. " " .and. i .ne. dStr%length .and. dStr%data(i:i) .ne. '' ) then
																foundChar = .true.
																iChar = iChar + 1
																dArr(j)%data(iChar:iChar) = dStr%data(i:i)
															elseif(dStr%data(i:i) .eq. " " .and. foundChar ) then
																dArr(j)%length = iChar
																foundChar = .false.
																iChar = 0
																j = j + 1
															elseif(i == dStr%length .and. dStr%data(i:i) .ne. " ") then
																iChar = iChar + 1
																dArr(j)%data(iChar:iChar) = dStr%data(i:i)
																dArr(j)%length = iChar
															endif
														enddo

														iCount = j
			!do i = 1, j
			!	dTmp = new(" dArr(") + intToStr(i) + new("): |") + dArr(i) + new("|")
			!	print *, dTmp%data(1:dTmp%length)
			!enddo	
														return
													end subroutine createStrArr

													subroutine pall(char,iNum)
														implicit none

            !Args
														character char
														integer iNum        

            !Locals
														type(dstring) :: dPall

														dPall = all(char,iNum)

														print *, dPall%data(1:dPall%length)
														return
													end subroutine pall

													function all(char,iNum)
														implicit none

        	!Function type
														type(dString) all

        	!Args
														character char
														integer iNum

        	!Locals
														integer i

														all%length=0 		        	
														!all%data(i:iNum) = ''		 		

														do i = 1, iNum
															all = all +	new(char)	
														enddo

														return	
													end function all

													function fit(dStr, dIn)
														implicit none

            ! Function type
														type(dString) :: fit

            ! Args
														type(dString) :: dStr, dIn

            ! Locals
														integer i, remains, iExtra
														real ratio
														type(dString) :: dBefore, dAfter

														remains = dStr%length - dIn%length
														ratio = real(remains)/2.

														iExtra = closeInt(ratio)

														dBefore = new(dStr%data(1:iExtra))
														dAfter = new(dStr%data(dStr%length-iExtra:dStr%length))

														fit = dBefore + dIn + dAfter
														return
													end function fit

													function trimUp(str)
														implicit none

          ! Type
														type(dstring) trimUp

          ! Args
														character (len=*) str

          ! Locals
														integer i
														integer myLen
														integer start, end
														integer diff

														myLen = len(str)

														do i = 1, myLen

															if(str(i:i) /= ' ') then
																start = i
																exit
															endif

														enddo

														do i = myLen, 1, -1

															if((str(i:i) /= ' ') .and. (str(i:i) /= '') ) then
																end = i
																exit
															endif 

														enddo

													diff = end - start

													trimUp%data(1:diff) = str(start:end)
													trimUp%length = diff

													return

													end function trimUp

													subroutine askteriskToPrime(string)
														implicit none

          ! Args
														character (len=*) string

          ! Locals
														integer i
														integer myLen

														myLen = strLen(string)

														do i = 1, myLen

															if(string(i:i) == "*") string(i:i) = "'"

														enddo

														return

													end subroutine askteriskToPrime


													subroutine replace_dstr(mainDstr, subDstr1, subDstr2)

														implicit none

														type(dstring) mainDstr, subDstr1, subDstr2

														call replace_str(mainDstr%data,subDstr1%data(1:subDstr1%length),subDstr2%data(1:subDstr1%length))
														mainDstr%length = strLen(mainDstr%data)

														return

													end subroutine replace_dstr

													subroutine replace_str(mainStr,subStr1,subStr2)
														implicit none

          ! Args
														character (len=*) mainStr
														character (len=*) subStr1,subStr2

          ! Locals
														integer iEndOfSearch
														integer iLenSub1, iLenSub2
														integer iLongestStr
														integer i
														integer occursAt(100)
														integer numOccurs

          ! Find the lengths of both strings
														iLenSub1 = len(subStr1)
														iLenSub2 = len(subStr2)

														iLongestStr = max(iLenSub1,iLenSub2)


          ! Find where the last place you should look for the occurances is
														iEndOfSearch = strLen(mainStr) - iLongestStr

														numOccurs = 0

          ! First, find all occurances
														do i = 1, iEndOfSearch

             !if(s_eqi_str(mainStr(i+iLenSub1-1),subStr1)) then ! found an occurance
             !   numOccurs = numOccurs + 1
             !   occursAt(numOccurs) = i
             !endif

														enddo



          ! Next, replace all occurances if there are any

														if (numOccurs > 0) then

															do i = 1, numOccurs

                ! Finish this later if you want....
                !mainStr(occursAt(i):iLenSub2

															enddo

														endif

													end subroutine replace_str

													function inStr(cMain,cSub,iPos)
														implicit none

			!Function type
														logical inStr

			!Args
														character (len=*) cMain,cSub
														integer, optional, intent(out) :: iPos

			!Locals 
														integer i
														type(dstring) dMain, dSub

														dMain = new(trim(cMain))
														dSub  = new(trim(cSub))

														inStr = .false.
														do i = 1, dMain%length   				
															if(i + (dSub%length - 1) > dMain%length) then 
																exit
															else
																if(s_eqi_str(dMain%data(i:i+(dSub%length - 1)),&
																	dSub%data(1:dSub%length))) then
																	inStr = .true.
																	if(present(iPos)) iPos = i
																	exit
																endif
															endif   				
														enddo

														return								
													end function inStr

													function indStr(dMain,dSub,iPos)
														implicit none

			!Function type
														logical indStr

			!Args
														type(dstring) :: dMain,dSub
														integer, optional :: iPos

														if(present(iPos)) then
															indStr = inStr(dMain%data(1:dMain%length),dSub%data(1:dMain%length),iPos)		
														else
															indStr = inStr(dMain%data(1:dMain%length),dSub%data(1:dMain%length))		
														endif	

														return								
													end function indStr



													function isNull(string)

          ! Function type
														logical isNull

          !Argument variables
														character (len=*) string

          !Local variables
														integer sLen
														logical myTest

														sLen = len(string)

														myTest = .true.

														do i = 1, sLen

															if(string(i:i) /= '')then

																myTest = .false.
																exit

															endif

														enddo

          !Return the logical value
														isNull = myTest
														return
													end function isNull

													subroutine trimPeriod(string)
														implicit none

             ! args
														type(dstring) string

             ! locals
														integer iStrLen, i


														iStrLen = string%length

														do i = iStrLen, 1, -1

															if (string%data(i:i) .ne. '' .or. string%data(i:i) .ne. ' ') then

																if ( string%data(i:i) == '.') then
																	string%data(i:i) = ''
																	string%length = iStrLen - 1
																	exit

																endif

															endif

														enddo

														return

													end subroutine trimPeriod

													subroutine trimAsterisk(string)
														implicit none

             ! args
														character (len=*) string

             ! locals
														integer iStrLen, i


														iStrLen = len(trim(string))

														do i = iStrLen, 1, -1

															if (string(i:i) .ne. '' .or. string(i:i) .ne. ' ') then

																if ( string(i:i) == '*') then
                      !trimAsterisk = .true. ! found an asterisk
																	string(i:i) = ''

																	exit

																endif

															endif

														enddo

														return

													end subroutine trimAsterisk

													subroutine getFileName(file,Name,Ext)
														implicit none

             ! Args
														character (len=*) file
														type(dstring) Name, Ext

             ! Locals
														integer myLen
														integer iDot


														myLen = strLen(file)
														iDot = index(file, ".")

														if(iDot > 0) then
															Name%data = file(1:iDot-1)
															Name%length = iDot - 1
															Ext%data = file(iDot+1:myLen)
															Ext%length = myLen - iDot
														else
															Name%data = file
															Name%length = myLen
															Ext%data = ''
															Ext%length = 0
														endif

														return
													end subroutine getFileName

													subroutine moveOver(string)
														implicit none

             !Args
														character ( len=*) string

             !Locals
														integer iStrLen, i, j

														character (len=125) buffer

														iStrLen = len(string)


														buffer(i:125) = ''

														do i = 1, iStrLen

															if ( string(i:i) .ne. '' .and. string(i:i) .ne. ' ') then
																exit

															endif

														enddo

														do j = 1, iStrLen

															if (i .gt. iStrLen) exit
															buffer(j:j) = string(i:i)
															i = i + 1

														enddo

														string(1:iStrLen) = ''
														string(1:iStrLen) = buffer(1:iStrLen)

													end subroutine moveOver

													function Es_T(temp)

          ! function type
														real Es_T

          !args
														integer :: temp ! in degrees K

														Es_T = 78.54 + (-0.3815)*(temp - 298.0)

														return

													end function Es_T

													subroutine suggestFileName_dstr(oldFile,newFile)

														implicit none

          ! Args
														type(dstring) oldFile, newFile

														call suggestFileName_str(oldFile%data(1:oldFile%length),newFile%data)
														newFile%length = strLen(newFile%data)

													end subroutine suggestFileName_dstr

													subroutine suggestFileName_str(oldFile, newFile)
														implicit none
          !Argument variables
														character (len=*) oldFile
														character (len=*) newFile


          !Local variables
														integer posDot
														integer posUnderS
														integer fileNumber
														integer numDigits
														integer i
														integer num
														integer myPos
														integer tempFileLen
														type(dstring) strNum
														character (len=10) fileExt
														character (len=100) fileName
														character charDot
														character (len=100) tempFile
														logical keepGoing

														tempFile = oldFile
														newFile = oldFile


														num = 0

														inquire(file=trim(tempFile), exist=keepGoing)
														i = 0
														do while(keepGoing)

															i = i + 1
															tempFileLen = strLen(tempFile)

															fileExt = ''

															posDot = index(tempFile, '.')
															posUnderS = index(tempFile, '__')

															if(posDot > 0) then 
																charDot = '.' 
															else

																charDot = ''

															endif

															if(posDot == 0 ) then ! There is no "."

																if(posUnderS > 0) then

																	num = i
																	strNum = intToStr(num)
																	fileName = tempFile(1:posUnderS-1)

																else ! It is the first time a file has been suggested

																	strNum%data = '1'
																	strNum%length = 1
																	fileName = tempFile(1:tempFileLen)

																endif ! End search for underscore

															else ! there is a ".", and probably a file extension

																fileExt = tempFile(posDot+1:tempFileLen)

																if(posUnderS > 0) then

																	num = i
																	strNum = intToStr(num)
																	fileName = tempFile(1:posUnderS-1)

																else ! It is the first time a file has been suggested

																	strNum%data = '1'
																	strNum%length = 1
																	fileName = tempFile(1:posDot-1)

																endif ! End search for underscore

															endif ! End check for "."

															newFile = trim(fileName) // '__' // strNum%data(1:strNum%length) // charDot // fileExt

															inquire(file=trim(newFile), exist=keepGoing)

															if(keepGoing .eqv. .true.) tempFile = newFile

														enddo

														return
													end subroutine suggestFileName_str


													function isNeg(int)
														implicit none

          ! Function type
														logical isNeg

          ! Args
														integer int

          ! Local
														logical temp

														if(int < 0) then

															temp = .true.

														else
															temp = .false.

														endif

														isNeg = temp
														return

													end function isNeg

													function logicalToStr(myBool)
														implicit none

          ! Function type
														type(dstring) logicalToStr

          ! Args
														logical myBool

          ! Locals
														integer numDigits
														integer iTemp
														integer i
														character (len=20) tmpStr

														tmpStr = ''
														iTemp = getUnit()
														open(unit=iTemp,status='SCRATCH')
														write(iTemp,'(l10)') myBool
														rewind(iTemp)
														read(iTemp, *) tmpStr
														close(iTemp)

          !numDigits = strLen(tmpStr)
														logicalToStr = new(trim(tmpStr))
          !logicalToStr%length = numDigits
          !logicalToStr%data(1:numDigits) = tmpStr(1:numDigits)

														return

													end function 							

													function realTodStr(myReal,nDigits)
														implicit none

          ! Function type
														type(dstring) realTodStr

          ! Args
														real myReal
														integer nDigits

          ! Locals
														integer numDigits
														integer iTemp
														integer i
														character (len=20) tmpStr
														logical lInStr
														integer iInStr
														type(dstring) :: dTmp

														tmpStr = ''
														iTemp = getUnit()
														open(unit=iTemp,status='SCRATCH')
														write(iTemp,'(f6.3)') myReal
														rewind(iTemp)
														read(iTemp, *) tmpStr
														close(iTemp)

          ! Erase any extra zeros
														do i = 20, 1, -1
															if(.not. s_eqi_str(tmpStr(i:i), '')) then
																if(s_eqi_str(tmpStr(i:i), '0')) then
																	tmpStr(i:i) = ''
																else
																	exit
																endif
															endif
														enddo

          !numDigits = strLen(tmpStr)

														realTodStr%length = nDigits
														realTodStr%data(1:nDigits) = tmpStr(1:nDigits)
														dTmp = realTodStr															
														do i = nDigits, 1, -1
															if( dTmp%data(i:i) .eq. " " ) dTmp%data(i:i) = "0"
														enddo

														realTodStr = dTmp

														return

													end function realTodStr


													function real8TodStr(myReal)
														implicit none

          ! Function type
														type(dstring) real8TodStr

          ! Args
														real*8 :: myReal

          ! Locals
														integer numDigits
														integer iTemp
														integer i
														character (len=80) tmpStr   
														type(dstring) :: dTmp


														dTmp = new("'(D)'") 
														iTemp = getUnit()
														open(unit=iTemp,status='SCRATCH')
														write(unit=iTemp,fmt='(D17.10)') myReal
														rewind(iTemp)
														read(iTemp, *) tmpStr
														close(iTemp)

														real8TodStr = new(tmpStr)													

														return
													end function real8TodStr  

													function tailInt(str)
														implicit none

          ! Function type
														integer tailInt

          ! Args
														character(len=*) str

          ! Locals
														integer iTemp
														integer iStart
														character (len=80) tmpStr
														type(dstring) :: dTmp

														dTmp = new(str)

														iStart = dTmp%length - 5

														dTmp = new(dTmp%data(iStart:dTmp%length))

														iTemp = getUnit()
														open(unit=iTemp,status='SCRATCH')
														write(iTemp,*) dTmp%data(1:dTmp%length)
														rewind(iTemp)
														read(iTemp, *) tailInt
														close(iTemp)

														return

													end function tailInt

													function realToStr(myReal)
														implicit none

          ! Function type
														type(dstring) realToStr

          ! Args
														real myReal

          ! Locals
														integer numDigits
														integer iTemp
														integer i
														character (len=20) tmpStr

														tmpStr = ''
														iTemp = getUnit()
														open(unit=iTemp,status='SCRATCH')
														write(iTemp,*) myReal
														rewind(iTemp)
														read(iTemp, *) tmpStr
														close(iTemp)

          ! Erase any extra zeros
														do i = 20, 1, -1
															if(.not. s_eqi_str(tmpStr(i:i), '')) then
																if(s_eqi_str(tmpStr(i:i), '0')) then
																	tmpStr(i:i) = ''
																else
																	exit
																endif
															endif
														enddo

														numDigits = strLen(tmpStr)

														realToStr%length = numDigits
														realToStr%data(1:numDigits) = tmpStr(1:numDigits)

														return

													end function realToStr


													function intToStr(int)
														implicit none

          ! Function type
														type(dstring) intToStr

          ! Args
														integer int

          ! Locals
														integer numDigits
														integer iTemp
														character (len=20) tmpStr


														iTemp = getUnit()
														open(unit=iTemp,status='SCRATCH')
														write(iTemp,*) int
														rewind(iTemp)
														read(iTemp, *) tmpStr
														close(iTemp)

														numDigits = strLen(tmpStr)

														intToStr%length = numDigits
														intToStr%data(1:numDigits) = tmpStr(1:numDigits)

														return

													end function intToStr

													function intToStr_long(int)
														implicit none

          ! Function type
														type(dstring) intToStr_long

          ! Args
														integer*8 int

          ! Locals
														integer numDigits
														integer iTemp
														character (len=20) tmpStr


														iTemp = getUnit()
														open(unit=iTemp,status='SCRATCH')
														write(iTemp,*) int
														rewind(iTemp)
														read(iTemp, *) tmpStr
														close(iTemp)

														numDigits = strLen(tmpStr)

														intToStr_long%length = numDigits
														intToStr_long%data(1:numDigits) = tmpStr(1:numDigits)

														return

													end function intToStr_long


													function strToInt(str)
														implicit none

          ! Function type
														integer strToInt

          ! Args
														character (len=*) str

          ! Locals
														integer pos 
														integer iStrLen
														integer i,j
														integer iExp
														integer iIntLen
														integer iDecLen
														integer iTemp
														character (len=20) intPart
														character (len=20) decPart
														logical isNegative
														real rTemp, temp



          ! Find the length of the string
														iStrLen = strLen(str)

          ! Set up initial conditions
														iExp = iStrLen
														if(str(1:1) == "-") then 
															isNegative = .true.
															intPart = str(2:iStrLen)

															iStrLen = iStrLen - 1

														else
															isNegative = .false.
															intPart = str(1:iStrLen)

														endif

          ! loop through only the integer part and sum value
														do i = 1, iStrLen
															iExp = iExp - 1
															iTemp = intCharToInt(intPart(i:i))

															rTemp = (iTemp*10**(iExp)) + rTemp

														enddo



														if (isNegative) rTemp = -rTemp

														strToInt = int(rTemp)

														return

													end function strToInt

													function strToReal(str)

														implicit none

          ! Function type
														real strToReal

          ! Args
														character (len=*) str

          ! Locals
														integer pos 
														integer iStrLen
														integer i,j
														integer iExp
														integer iIntLen
														integer iDecLen
														integer iTemp
														character (len=20) intPart
														character (len=20) decPart
														logical isNegative
														real rTemp, temp


          ! Find the position of the decimal
														pos = index(trim(str), ".")

														if (pos > 0 ) then
             ! Find the length of the string
															iStrLen = strLen(str)

             ! Set up initial conditions
															if(str(1:1) == "-") then 
																isNegative = .true.
																intPart = str(2:pos - 1)


															else
																isNegative = .false.
																intPart = str(1:pos - 1)


															endif
															decPart = str(pos+1:iStrLen)

             !Find the true lengths
															iIntLen = strLen(intPart)
															iDecLen = strLen(decPart)

             ! More initial conditions
															rTemp = 0.
															iExp = iIntLen

             ! Loop through to sum the integer part
															do i = 1, iIntLen
																iExp = iExp - 1
																iTemp = intCharToInt(intPart(i:i))

																rTemp = (iTemp*10**(iExp)) + rTemp

															enddo


             ! Loop through to sum the decimal part
															do i = 1, iDecLen


																iTemp = intCharToInt(decPart(i:i))

																temp = 10**(i)

																temp = 1/temp
																rTemp = iTemp*temp + rTemp

															enddo

														else ! There is no '.' in the string, so you are just dealing with an integer

             ! Find the length of the string
															iStrLen = strLen(str)

             ! Set up initial conditions
															iExp = iStrLen
															if(str(1:1) == "-") then 
																isNegative = .true.
																intPart = str(2:iStrLen)

																iStrLen = iStrLen - 1

															else
																isNegative = .false.
																intPart = str(1:iStrLen)

															endif

             ! loop through only the integer part and sum value
															do i = 1, iStrLen
																iExp = iExp - 1
																iTemp = intCharToInt(intPart(i:i))

																rTemp = (iTemp*10**(iExp)) + rTemp

															enddo

														endif

														if (isNegative .eqv. .true.) then
                                 rTemp = -rTemp
                                endif

														strToReal = rTemp

														return

													end function strToReal


													function strLen(str)
														implicit none

          ! Function type
														integer strLen

          ! Args
														character (len=*) str

          ! Locals
														integer i
														integer iStrLen

														iStrLen = len(str)

														do i = iStrLen, 1, -1

															if((str(i:i) .ne. ' ') .or. (str(i:i) .ne. '')) exit

														enddo

														strLen = i
														return

													end function strLen

													function intToIntChar(int)
														implicit none

          ! Function type
														character intToIntChar

          !Argument variables
														integer int

          !Local variables
														character (len=1) iChar

														if (int == 1) then 
															iChar = '1'
														elseif(int == 2)then
															iChar = '2'
														elseif(int == 3)then
															iChar = '3'
														elseif(int == 4)then 
															iChar = '4'
														elseif(int == 5)then
															iChar = '5'
														elseif(int == 6)then
															iChar = '6'
														elseif(int == 7)then
															iChar = '7'
														elseif(int == 8)then
															iChar = '8'
														elseif(int == 9)then
															iChar = '9'
														elseif(int == 0)then
															iChar = '0'
														endif

														intToIntChar = iChar

														return

													end function intToIntChar

													function intCharToInt(char)

														implicit none

          ! Function type 
														integer intCharToInt

          !Argument variables
														character (len=1) char

          !Local variables
														integer iTemp

														if ( char == '1') then
															iTemp = 1
														elseif(char == '2')then
															iTemp = 2
														elseif(char == '3')then
															iTemp = 3
														elseif(char == '4')then 
															iTemp = 4
														elseif(char == '5')then
															iTemp = 5
														elseif(char == '6')then
															iTemp = 6
														elseif(char == '7')then
															iTemp = 7
														elseif(char == '8')then
															iTemp = 8
														elseif(char == '9')then
															iTemp = 9
														elseif(char == '0')then
															iTemp = 0
														endif

														intCharToInt = iTemp

														return

													end function intCharToInt

													function isInteger(dstr)
														implicit none

	      ! Function type 
														logical isInteger

          ! Argument variables
														type(dstring) dstr

														integer i


														isInteger = .true.

														do i = 1, dstr%length
															if(dstr%data(i:i) .eq. ".") then
																isInteger = .false.
																exit
															endif
														enddo 

          !write(*,'(a)') " dstr: ", dstr%data(1:dstr%length)

          !write(*,*) " isInteger = ", isInteger

														return
													end function isInteger



													function isNumeric(char)
														implicit none

          ! Function type 
														logical isNumeric

          ! Argument variables
														character ( len = 1 ) char

														isNumeric = .false.

														if ( char == '1' .or. &
															char == '2' .or. &
															char == '3' .or. &
															char == '4' .or. &
															char == '5' .or. &
															char == '6' .or. &
															char == '7' .or. &
															char == '8' .or. &
															char == '9' .or. &
															char == '0') then

															isNumeric = .true.
														endif

														return
													end function isNumeric

													function closeInt(myReal)
														implicit none

         	! Function type
														integer closeInt

         	! Args
														real myReal										

         	! Locals
														integer myInt
														real fromTop, fromBottom
														integer one

														if (myReal < 0.0 ) then
															one = -1
														else
															one = 1
														endif

														myInt = myReal


														fromTop = abs((myInt + one) - myReal)
														fromBottom = abs(myInt - myReal)

														if (fromTop < fromBottom) then
															closeInt = myInt + one
														else
															closeInt = myInt
														endif	

														return
													end function closeInt


													function lmod_int(i1,i2)
														implicit none

          ! Function type 
														logical :: lmod_int

          ! args
														integer :: i1, i2

														if(mod(i1,i2) .eq. 0) then
															lmod_int = .true.
														else
															lmod_int = .false.
														endif

														return
													end function lmod_int

													function lmod_real8(i1,i2)
														implicit none

          ! Function type 
														logical :: lmod_real8

          ! args
														real*8 :: i1, i2


														if(mod(i1,i2) .gt. 0.d0) then
															lmod_real8 = .false.
														else
															lmod_real8 = .true.
														endif

														return
													end function lmod_real8

													function lmod_real4(i1,i2)
														implicit none

          ! Function type 
														logical :: lmod_real4

          ! args
														real :: i1, i2


														if(mod(i1,i2) .gt. 0.d0) then
															lmod_real4 = .false.
														else
															lmod_real4 = .true.
														endif

														return
													end function lmod_real4

													function str2bool(string)
														implicit none

          ! Function type 
														logical str2bool

          ! args
														character (len=*) string

          !if ( s_eqi_str(string, 'Y') .or. &
          !     s_eqi_str(string, '1') .or. &
          !     s_eqi_str(string, 'T') .or. &
          !     s_eqi_str(string, '.true.')) then
														if ( s_eqi_str(string(1:1), 'Y') .or. &
															s_eqi_str(string(1:1), '1') .or. &
															s_eqi_str(string(1:1), 'T') .or. &
															s_eqi_str(string, '.true.')) then

															str2bool = .true.

														else
															str2bool = .false.

														endif

														return


													end function str2bool

													function s_eqi_dstr(dstr1,dstr2)

														implicit none

          ! Function type
														logical s_eqi_dstr

          ! Args
														type(dstring) dstr1, dstr2

														s_eqi_dstr = s_eqi_str(trim(dstr1%data(1:dstr1%length)), &
															trim(dstr2%data(1:dstr2%length)))

														return

													end function s_eqi_dstr

													function s_eqi_str ( strng1, strng2 )
          ! Function type 
														logical s_eqi_str

														character (len=*)  strng1, strng2

														integer i
														integer len1
														integer len2
														integer lenc
														character s1
														character s2
          !
														len1 = len ( strng1 )
														len2 = len ( strng2 )
														lenc = min ( len1, len2 )

														s_eqi_str = .false.

														do i = 1, lenc

															s1 = strng1(i:i)
															s2 = strng2(i:i)
															call c_cap ( s1 )
															call c_cap ( s2 )

															if ( s1 /= s2 ) then
																return
															end if

														end do

														do i = lenc + 1, len1
															if ( strng1(i:i) /= ' ' ) then
																return
															end if
														end do

														do i = lenc + 1, len2
															if ( strng2(i:i) /= ' ' ) then
																return
															end if
														end do

														s_eqi_str = .true.

														return
													end function s_eqi_str

													subroutine c_cap ( c )

          ! Argument variable
														character c

          ! Local variables
														integer itemp

														itemp = ichar ( c )

														if ( 97 <= itemp .and. itemp <= 122 ) then
															c = char ( itemp - 32 )
														end if

														return
													end subroutine c_cap

													function getUnit()
														implicit none
         ! Function type 
														integer getUnit

														integer i
														integer ios
														integer iunit
														logical lopen



														iunit = 0

														do i = 1, 99

															if ( i /= 5 .and. i /= 6 ) then

																inquire ( unit = i, opened = lopen, iostat = ios )

																if ( ios == 0 ) then
																	if ( .not. lopen ) then
																		getUnit = i
																		exit
																	end if
																end if

															end if

														end do

														return

													end function getUnit

													subroutine warn(sFuncName, sErr)

            ! Declare local variables
														character (len=*) sFuncName
														character  (len=*) sErr
														type(dString) :: dWarn, dTopBottom

														dWarn = new(' Warning[') + new(' ') + new(sFuncName) + new(' ->') + new(' ') + new(sErr) + new(' ]')
														dTopBottom = all("=", dWarn%length)
														write (*,'(a)') dTopBottom%data(1:dTopBottom%length)
														write (*,'(a)') dWarn%data(1:dWarn%length)
														write (*,'(a)') dTopBottom%data(1:dTopBottom%length)

													end subroutine warn


													subroutine raiseError(sFuncName, sErr)
														implicit none
            ! Declare local variables
														character (len=*) sFuncName
														character  (len=*) sErr


														type(dString) dWarn, dTop, dBottom, dInsert

														dWarn = new('| FATAL ERROR : [') + new(' ') + new(sFuncName) + new(' ->') + new(' ') + new(sErr) + new(' ] |')
														dTop = all("=", dWarn%length)
														dInsert = new("Program Ended Abruptly")
														dBottom = fit(dTop, dInsert)

														if(dBottom%length .gt. dTop%length) dTop = all("=", dBottom%length)
														write (*,*) ""
														write (*,'(a)') dTop%data(1:dTop%length)
														write (*,'(a)') dWarn%data(1:dWarn%length)
														write (*,'(a)') dBottom%data(1:dBottom%length)
														stop ' *** UNNATURAL PROGRAM STOP ***'
            !stop            dBottom%data(1:dBottom%length)


													end subroutine raiseError

												end module commonFunctions


