function Compile_CWT

MEX_OK = 0;

% for file = {'fwt','ifwt','fwt2','ifwt2', 'fwt2_CWT', 'afwt2_CWT',  ...
%         'ifwt2_CWT', 'afwt','afwt2'}
%     
%     file = char(file);
%     if exist(file)~=3,
%         MEX_OK = 0;
%         break;
%     end
% end

if ~MEX_OK
    Friend = computer;
    isPC = 0;
    if strcmp(Friend(1:2),'PC')
        isPC = 1;
    end
    if isPC
        mex -c fwt_level.c
        mex fwt.c fwt_level.obj -output fwt
        mex ifwt.c fwt_level.obj -output ifwt
        mex fwt2.c fwt_level.obj -output fwt2
        mex ifwt2.c fwt_level.obj -output ifwt2
        
        mex -c afwt_level.c
        mex afwt.c afwt_level.obj -output afwt
        mex afwt2.c afwt_level.obj -output afwt2
        
        mex fwt2_CWT.c fwt_level.obj -output fwt2_CWT
        mex ifwt2_CWT.c fwt_level.obj -output ifwt2_CWT
        mex afwt2_CWT.c afwt_level.obj -output afwt2_CWT
        
        % mex aifwt_mod.c fwt_level.obj -output aifwt_mod
        % mex aifwt2_mod.c fwt_level.obj -output aifwt2_mod
    else       
        mex -v GCC='/usr/bin/gcc-4.9' -c fwt_level.c
        mex -v GCC='/usr/bin/gcc-4.9' fwt.c fwt_level.o -output fwt
        mex -v GCC='/usr/bin/gcc-4.9' ifwt.c fwt_level.o -output ifwt
        mex -v GCC='/usr/bin/gcc-4.9' fwt2.c fwt_level.o -output fwt2
        mex -v GCC='/usr/bin/gcc-4.9' ifwt2.c fwt_level.o -output ifwt2

        mex -v GCC='/usr/bin/gcc-4.9' fwt2_CWT.c fwt_level.o -output fwt2_CWT
        mex -v GCC='/usr/bin/gcc-4.9' ifwt2_CWT.c fwt_level.o -output ifwt2_CWT

        mex -v GCC='/usr/bin/gcc-4.9' -c afwt_level.c
        mex -v GCC='/usr/bin/gcc-4.9' afwt.c afwt_level.o -output afwt
        mex -v GCC='/usr/bin/gcc-4.9' afwt2.c afwt_level.o -output afwt2
        mex -v GCC='/usr/bin/gcc-4.9' afwt2_CWT.c afwt_level.o -output afwt2_CWT
        
        % mex aifwt.c fwt_level.o -output aifwt_mod
        % mex aifwt2.c fwt_level.o -output aifwt2_mod
    end
end
disp('Done compiling!')