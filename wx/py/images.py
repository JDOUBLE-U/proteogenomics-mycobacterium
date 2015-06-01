"""Support for icons."""

__author__ = "Patrick K. O'Brien <pobrien@orbtech.com> / David Mashburn <david.n.mashburn@gmail.com>"

import wx
from wx.lib.six import BytesIO

def getPyIcon(shellName='PyCrust'):
    icon = wx.Icon()
    icon.CopyFromBitmap(getPyBitmap(shellName))
    return icon

def getPyBitmap(shellName='PyCrust'):
    return wx.Bitmap(getPyImage(shellName))

def getPyImage(shellName='PyCrust'):
    stream = BytesIO(getPyData(shellName))
    return wx.Image(stream)

def getPyData(shellName='PyCrust'):
    if shellName=='PyCrust':     
        return \
b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00 \x00\x00\x00 \x08\x06\x00\
\x00\x00szz\xf4\x00\x00\x00\x04sBIT\x08\x08\x08\x08|\x08d\x88\x00\x00\x04\
\x95IDATx\x9c\xed\x97?lSG\x1c\xc7?\x97\x98\xd8Q\xa3\xdeY\xa2j\x06\xa4\xf7"QJ\
\xbb<3@\x01\xa9\xc2\x0c\xa8!\x1d\x1c6\xcaB\xa8D[uI2\xf4\x8f\xe8\x103\xb4\xa2\
,5\x0b\x03\x032C\xab\xc0\x92dh:t\xc0)*E\xcd@<Q\x01Rl\t\xd4D\xaa\xe4{R\xd0{&\
\xa5\xd7\xe1\xfc\xec\xe7\xfciR\x08e\xe9O\xb2\xee|\xfe\xbd\xfb}~\xdf\xdf\xbd\
\xbb\xb3PJ\xf1"\xad\xe3\x85F\xff\x1f\xe0y\x03h\xad\xcdA\xc7~\xb4\xd6f-\x9f\
\xc4\xf3\x0c>Y\x1c#\x97\xddCUk\xf4B\x8d3\x9f\x8d\x9a\x9bU%\xe2~b\xab\xdf\x82\
\x83N+\xd3\xe92\\\x1f\xcf\x93\xdd\x9f\xa1\xaa5\x95\xf9\n\xe7\xf3y\xe2\x10[V\
\x82H\xe6\xd3G\x1dN\xf7\xc3\xa7\xc7a\xc0\x83\xc3\xc7\xf3\xcc\xcc\xcd\xe3(\
\x85\xdb\xe7\xf2\xc9\xe8X\x1b\xe43+\x10\xd5\xb6\x94\x87Z\xe8\x90NU\x91I@\x00\
\x06\xbe\x18\xb7J\x98\xca$`\x98\xb9]&{,\x8fRV\x85\xa7V@k\x9bq)o\x83+\t\xe9T\
\xd5f\x95\x02\x91\xb4~_\r\xd9\xb6\xbaP\x03\x04n\x9f\xcbDa\xb8\t\xfe\xaf\x17a\
<\xe3\xc8\x94lo\x9b\xd6\xa8\xf4\x80\x07\xb7o\xcd\xe0\x0c\x0e\xa2R\x8a\xb4\
\x93n\xbal\x1a`e\xe0U\xc1\xd6\xb0\xb8\n\x99\x91"\x93\xaf\xba\xe4\x0ed\xda|6,\
\x81\xd6\xda\x9c|\xab]\xea\xcd\x04\x8f\x9b\t\xad\nz\xa1\x02\x80\xdb\xe7R\x1a\
\xcf\xa3\xb56\xeb\x02D5\x9e\xf8\xdc\xe1T\xff\xd3\x05\x8e\x82\x83U\xe1Z\xb1\
\x18\x9b\xbf\x06\xacQ\x82H\xea\x01/Z@Ut\x08R\xb4$}\x16\xd3\x81A&%\xde\xee\
\xbev\x80x\xe0]{\xb2\x1cR\xa5\xe6C*\xb5\xf1\xc4Q\xa6"e\xfbQ\x1b\x8dE\xe6\x87\
>\xaa[Q\xadi\x0b\xb0r\x8f\x9e.\xc3t\xb9\xc4]\xaf5\xf6\xfe\xdb\xddt&\x02\xfa\
\x9c\xf5\x01\xe2A\xa2\xbeX\x01>]ntR\x12\xe3[\x00\x01\x98\x89\x11[_\xed\xafn\
\xab\x81U\xa0\xe7I7\x00\x97o\x04\xcd\x89\x06<;\xe9\x80\x07]i\x97\xc17\x1f\
\xd2\xd3\x91`\xe9\xaf?\x01p^Y\x06Z\n\xfau8?a\xfb]i\x97\xec\xa1\x8c\x05(|\xd8\
N\xba\xb3\xab\x87\xfb\x8f\x97\xd8\xd9\xd5\x03\xc0\xfd\xc7K\xec\xd8\xd6\xdd\
\xfc\xfd\xc1r\xd0\xf4\x01\xda~\x03H\xf4\x04\xd4j :\xb75\xc7\xae\xfd\xbcLW\
\xda\xa5\xf0M\x1e\t\xcc\xcc\xcdq\xa9P@\x8c\xf5fL\xdaHF\x16g\x9a\x19\xad\xcc\
\xee\xcb\xa3\n\xad\xa1\xda\xf1\x08\xef\xe5\x97x\xf8\xc8f\xf8\xc7\x93:\xdb;\
\x93M\xc8\x08j\xc7\xb6n\x1e,\x07m`\x97o\x04|;>\xd1T\xc4\x17\x8a\x13\xb9\xc3\
\x88\x01\x0fs\xa4\x9cc}\xf3A\x190\x82\x1f\xddR{-\x1bV\xfc\xd8f\xba\xbd3\xd9\
\x06\x15\x07\xbb\xf8\xd3\x12\xdf]-"\x93\xb2\xb1C*\xde\xcd\x1d\xde\xccN(\xc1\
\xae\x17"\xd0#+<j\x17m{\xcd\x9bj\x00.\xaf\xf0Xb\xb8\xdfA\xa6\x14\x18\x03\x06\
\xb4o\xcf\x8d\xc4\xbervc\x86M\xdaz\x80\x00\x95T\x19?\xd0 @&%~c\xbc\xe3W\xaf\
\xb4e\x00\xffh\xc6@\xbd\x11\xbc\xde\x1a\xfe\xef.\xa5\xa2q4\n0\x81\xad\xe9\
\xae7<\x12\xaf\xf5\xc2hy\xaa\xe97\x9cS\\\x98\xb2\x0e\x03\xb1\xcdhW\xdaC\x1a\
\xa0\xa2\xa0\x0e"\x14`\xb0Y\x85\x1b\x1f\x12\xaa7\x03)\xd9\x84\xa8\xccW\xb8{\
\xa7L\xe2\xde\x02\x94\xc6Gp_\xcf\x80\x90\x98\xd0g\xf4\xac\x82Pc\x82\x1a\xd5\
\x10\x08}&\xa7J\xcc\xde.1]m\x80\xf6+\xee\xfd\xae\x9bo\xc4\xf0;\x80\xef\x90\
\x0e\x04\x06`Q!\x02\x05\xc2 \xb5\xc2\x95\x15d\xb4C&[\xf7\xd2\x04\x80\xbb\xdb\
\x9e\xd1\x8e\x02\x90\xd8\xd4$ I\x87\x80\xf1\xf1\xdc!4\xc3\x88\x94}\xd8,TH\
\xbb.5m\xf0C\x9f3\x1f\r\x01\x96.\x82\x1a9\xe9Q\xb8\xd2\xf8\xf25\x0c\xbe\xe7#\
\x92\x12\x1d[\x03\t\x00E\xf4\xa6\t\xaaZ7`$\x18\x90\xf8\xf8\x80JK\x94\xa1\x01\
\x07\xb8\x0e~X\xc3\xed\x16\xf8)\xf8~j\x12B\rI\x89_\xf7!0 \x04\xf9Q\xc0\x18\
\x0c\xd1i\xea\x13\xb7\x04\xc0\x89\x93C\xabj\xb6\xf7@\x96\xd9_J|0:\x86R\n\xb7\
\xd7@\xaa%\x9d\xa3$\xba.\x90RA]\xe3\x87\x1a\x89\xdd\xefeR\xc2\x1a\'\xa8\x1f\
\x82\x0e-@m\xd1\xde\x076\xbc\x15\x97~(\x9a\x89b\x9e\xd9[s\xab!\xf7g\xd6\x1c\
\x8f\xdb\xbel\x8e\xa1S\xc7\xda\xc6\xe6\xee\xccs\xe9\xdcYnV\x95\xd8\xf2?&q+\
\x9c\x1b1\xf3\xbf\xcd3{\xfdJ\xdb\xf8\xde\xfd\x19.\\\xad\x08\x80\xbf\x01\xd1\
\x86\xfa\x8b\xc7\xc0\xc8\xb7\x00\x00\x00\x00IEND\xaeB`\x82'
    elif shellName=='PySlices':
        return \
b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00 \x00\x00\x00 \x08\x06\x00\x00\
\x00szz\xf4\x00\x00\x00\x01sRGB\x00\xae\xce\x1c\xe9\x00\x00\x00\x06bKGD\x00\xff\
\x00\xff\x00\xff\xa0\xbd\xa7\x93\x00\x00\x00\tpHYs\x00\x00\x0b\x13\x00\x00\x0b\
\x13\x01\x00\x9a\x9c\x18\x00\x00\x00\x07tIME\x07\xd8\n\x16\x03#\x0eV,\xc0?\x00\
\x00\x06\xb8IDATX\xc3\xbd\x97\x7flS\xd7\x15\xc7?\xd7v\xe2\x17\x12x\xcf\x04X(\xb4\
vT\n\x94L}\xde\xd6\x1f\x012\x91Pe\xd0LZ\x82\x80\xf2\xa3\x1a\t\xd2\xc6\xd0\xa4)\
\xc9\xb4\x9f\x19\x1bN\xd6N\x1bRGZ\r\xb6vR\x15\xa6u\xb0I[b\x952\xf6C\xc2\xd01\xca\
"Z\xa7\x82\xb2-lq\x10\xac\x89`\xf2{4\x89\xed8\xf1\xdd\x1f\xd7v\xe2\xfc\xa2H\xed\
\x8ed\xbd\xfb\xae\xcf\xbd\xe7{\xce=\xdfs\xcf\x13\xcc,K\x97\x15\xf3\x9f#_\x04\xa3\
\xd0Ar\\b\xc5\x04\xef\xf4\xbb\xa8*\x1b\xe7\xd0\xeb\xc5\xa9\x84X\xc4\xdem\x15\x94\
,6p9\x1d\xa2x\xa1!*w>\x0b\xc0\xfaG\x1f>\xfa\xe7\xf3=G\xe2\xf1\xd1w\xb9\x8b\x88\x99\
&\x1f{\x90\xe1\xf2\x95D?\xf7\x04K\x8cB\xf2\x1c\x0ep:@\x02\x89\xa4@\xcb\x93\x9c\xbc\
\x94\xc7\xf9\xde"\xecQCn\xaaX!?\xfeP\tK\x97\xadHY\xc3Iy<x>\x15\xfc\xe3E\xc7\x9d\xf7\
G6\x00\x17\xee\x15\x80o\xdf\x93\\\xae-\xc7Y\xe2A\x13\x0e@\x82\xc3\x91\xab4:\x06\xe3\
)\xc1\xa0\r{\x8f\x14\x12\xb5\x87\xd8\xb3}\xe3\xf8W\xf7\xd5\x8e\xcf\x9b\xa7\xc5\xae\
\xfc\xfbfd\xcb\xe7\xdbL`3\xf0\x87{\x01\xb0rc\x19\xe1gw\x93\xafi8\xa7\x1a\x9e,\xef\x8f\
\xc0_\xae@b\x0c\xd6\xafv\xf2\xf4\x0b\x05D\xed!\xaew\xbf\xc2B\xbd\x88\xdb##\xf8\xfc\r\
\x00\x07\x81\xb6\x0f\x04\xa0r\r\xe3\xa6\x8f\xe4\x9e\x8d\xe4;\x1c\x08!f\x07p\xee2|\xaa\
\x14\x9c.\x88\xc7\xc1\x99\x0f\xd7n\xeb\x1c\xbb\xe8\xa5\xe5+;\xa8./\xa3\xdf\xb2x\xb7\
\xe7_\xd4\xd4\x7f?\x04TM\xddc\x9a\x7f+\xefCT\x9bs\x1b\xbfz\x03~z\x1a\x8a\xe7\x83p\x80\
\xe6\x06C\x87\xf9\x05\xf0\x89\xfbm\xda\xb7\xbd\xc3gv~\x87`\xe8m\xbc\x86\xc1\x1a\xf3A\
\xc2\xa7\x0fW\xd6\x98\xc8\xbb\x02\xb0c\xe0\xce\x03)\xd5o&y\xeb\x1al}\x0c\x96\x190\xaf\
`f\x9d\x1a\x13\xea\x1aZ\t\xbd\x19\xc6k\x18\x18%\x1e\xbe\x11\x080\x15\x84s\xf2\xcb\xaa\
\xfb\x90;\xd6#|\x1fC\xb8\\\n@&\n\x89$<\xfa5\xb8u\x07\x9e\xa9\x00-\x0f\x16,\x98\xfdx\x9e\
\xf4\xc3\xf5\x018p$DU\xe5\x13\x98\xbe\xa5\x88\x02\x8d5\xa5&\xb7\xfa\xce\x06z\x07i\x05pM^T\
\xb1\x1a\xca\x1e\x80\xfc\xbc)Q\x19\x81\x17O\xc2\x1bm\x8a\x8a.\xd7\xec\x9e\x03\xc8\xb8\xca\
\xae\xe7v\xa6\xf3\xaa\xae\x19\x19\xe9\xe2\xec\x990\xe1\x7f\xf4s\xaag\x86#x\xe6\xd3\xc8\xe5\
\xc50\xcf\x9d\xceN1A\xbdC\x9d\xb0c-\x8c\xa7\xd4y\x17hw).\x1a\x087\xfc\xee-\xb0\x13j\xee\x97\
\x9d!\xaa*L\x1a\xf7\xd7\xd2\xd9\xde\x88"\xf7\xa4\x08\xdc\xf8/|\xa1Z%\x95\x10*\xfc\xf6\x1d\
\xa8\x0c@(\x00yy\xca\xb8\xd39\xb7\xf1_\x9c\x81\x01\x1b\x12\xa3p`+l]\x07-\x1dP$,\x84t\xe0\
5t(\xf7g\xf5\xb3\x00\x16\xce\x87"\rR)H\x8e\xc3\xf3A\x05"\xf8Me\xd8\xe5\x9an\xfc\xf2ux\xf5\
\x1c\xfc\xfe\x12|{+\xdc\xbf\x18\xf6TM\x07\xf5\x83\x06\xf07u\xd0u\xdc\xc7\xf2%~"\xef\xf5\
\xe7\x1e\xc1S\x9fD.1\x94\x97c)\xd8\xdc\x06\x8f<\x00\r\x1b\xc0(\x84\x94T\x002\xf2\xbd\xe3\
\xe0o\x86A\x0b\xea\xab\xe0\xe4\x01\xd8Q\x01\xebV\xcd\x9e\x135&X\x03\x11UjK}\x84N\x04\x00\
\xa4\x0b@\xd7`\xd7Z\xa5\xfc\xd2i\x08~]\x8d\xc7Rp\xcb\x86\xfd?\x87M~\xe5qO\x04\xc2\x87\xa1\
m\x17\x1fHd|"\n-\x1d\x1d\xd4o\xa9\x03\xc0\xb2\xa2\xea\x08\xb6\xaf\xa3\xfb\xc4_\xe1[\xdb\x95\
W\xa1\x80Z\xb0\xed\xc7\xb0\xa1\x0c~{A\x19\xfc\xb0\xc4\x8aIt\xb7\x8e\xb9\xba4\x9d\x03)J\x1f\
\xf1*\xe3Ue\xf0\xfck\xf0\xda\xa5\t\xa3\xdf}\xfa\xde<\x15\x9a\x1ag\x9e\x99\xb9,\xa5\xe36F\
\x81A\x7f\xd4R\x00.\xf62\x1cKR\x0cp\xe6\xca\x84bK\xc7\xc48\xc3\xdb\xb9"1\xd9Hf,\xa6\xd05\
\xcb\x7fMG\xdaV\xf62\x92\x9dM\xaa\x96[\xf6\xf4g\x7f\x0c\x0c\r\x8a\xc6\x0b\xd8\xf6\x93\x185\
\xe6\xc4F\x99q\x8d\t\xf9\x1e\x1fuknR\xe4p1\x94\x1a\x03\xc0\xbb8\t\xa8\xf5\xa0j\xc2\xa1N5\
\xce\xf7\xf8\xa8\xdc\xe0W\x00\xda\xbf\x94\x8btE~\x11\xd7F\x87X\x91_\x04\xc0\xb5\xd1!\x96\
\xe7M\x94\xbe\x1b\xc9XV\x07\xc8\xf9\x0f\xc0U\x14#\x1a\x05\xe1\x9c(\xa9\xbf9\x9f$\xdf\xe3\
\xa3\xfdp\x00\x1d8\x1b\x0e\xf3r{;\xe2`\x89_z\xa4N\xd3\xe0\xd9Y\xbd{\xee)\x03\xcb\x82~\xc70\
\xe6\x82Bn\x0e+\x0fo\x8f\'X\xe4tgAf@-\xcf+\xe0F2\x96\x03\xec\x957b\xbcz\xa23\x1b\x11[\x18\
\xec\xae\xadB\xd4\x98\xc8\xea\x9e\xda9\xd2\xcb\x06C\x82\x14\xfc\xc9\x17\xca=\xcb\xb4t|Yy\xba\
\xc8\xe9\xce\x015\x19\xd8\xd1sC\xfc\xea\xd7\x1d\xe8n=\xdd\x85\x18|\xb6\xb6*\xf72\x9aYtP\
\xf9B\x06h\xf5\x14\x8d\xe8\xd1t\x98\xcd`\x1a`r\x8a\xc6\x10\x8d\x9b\xbd\xe8\x9a\x91\xbe\xe7\
\xc1\xb2U5t=\xdeS\xf9\xa1q|6\x80\x00\x11\xad\x07;f\x81\x00\xdd\xadcgJ\xf1\xdf\xcc\x10\xff\
\x17\x91\x12\x12i\xe3\x899:\xa2\x8fL2\x9d\x8d\x00\x19Sg\xba\xf2a\x13\xd7C%\xd0\xdc\x13\xcc\
\xea5\xd6\x1a\xbc\x10\xb4\xb2L\xc8\xf6\x8a\x1e\x13]\x02\x11\x03\x12 \xe2"\xfd\xa1\x00z\\\xbf\
\xab}\xa3\xc4\x0f\x9a\x9e\x05\x11\xe9\x8b\xf0\xcf\xab=\xb8z\x07 t\xa2\t\xdf*?\x08\x1d\x19\
\xb7in5 n!cQ\xfa\xe3@\xdc\xa6+\x18\xa2\xfb\xed\x10\xa7\xd27i\xe3f\x83\xde\xf7\xac,#\x1a7\x01\
\xb6\x17OL\xa8Nc\xd0@\xc4\x0c\x10\x12\xdd2\xf0\xe9\x11\xf4L\x85t\x8b\xdc~\xc0\xb7Z5\x08^#\x9d\
\xf5\xc8\xf4S\xc7\x13\x07\xa4\x8d\xe9k\xc0\xa2\x11\xa1\xa9\xc5r \x82\xc7\xe7#jI\xec\xb8M\xcb\
\xfe\x06\xa0?\x87\xa6M\xf5&\xed\xc7\xd2/?\x82\xba]6\xc2\xadcM\xca\x01\x97bd\x86i\x82~\xcbJ\x83\
\xd1A\x82\x8e\x8d\r\x18\x1e\x1dCf;)\xf0y\xb1\xe3Q|\x05\x02[\x83\xd7\x83]\x10\xb7\xc0\xadc\'l\x88\
\xa9\x8e6\xd0\xac\x12P\xa2\xda4\xb2\xf9?\t\xc0\xee\xfa\x86\xe9\xdf\x87k+\xe9\xbe\x10b_\xf3A\x0c\
\xc3\xc0W"A\x9b\x08\x9d\xd7\xd0\xb1\x12\x02]7 aa\xc7-\xf4t\xbd\xd7\xdd:\xcc\xd07\xdaq\xb0\xe2\
\n@t0:\xfb\xc7\xe9d\t_\xe8\x92\x9d\x1d\x01\xba\xdf\x0cO\x07Y\xee\x9fq~\xb2<^YK\xc3\xde-\xb9{^\
\xed\xe3\xe5\x1f\xb6r\xaa\x07!>J\xe6u\xfc, \xfb\xfe\xdeG\xf7\x99c\xd3\x80\xb7\xbe\x14\x16\x00\
\xff\x03\x07\x06\x80\xbbd\xd4\xb0\x14\x00\x00\x00\x00IEND\xaeB`\x82'
    elif shellName=='SymPySlices':
        return \
b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00 \x00\x00\x00 \x08\x06\x00\
\x00\x00szz\xf4\x00\x00\x00\x01sRGB\x00\xae\xce\x1c\xe9\x00\x00\x00\x06bKGD\
\x00\xff\x00\xff\x00\xff\xa0\xbd\xa7\x93\x00\x00\x00\tpHYs\x00\x00\x0b\x13\x00\
\x00\x0b\x13\x01\x00\x9a\x9c\x18\x00\x00\x00\x07tIME\x07\xd9\n\x10\x10\x0c\
\x119od\xf7\x00\x00\x07\xa2IDATX\xc3\xb5\x97{PT\xd7\x1d\xc7?wY\x96\xbb\xb0\xcb\
\xbd\x08\x08\xa8\x19V\x05\xe2\xab\xd95\xc1\x88\xa2\x03\x98\xf1\x11\'-\xda\x98h\
\x92\xa9b\x1fN:mG\xb0\xb13\xd5v\x80\xc4\xa9m\x98\x04\x8cQcfj \xcdL\x9c\xfc\
\x03$\xea\x18k\'\xc4\x91\xa8\xa5Iv\x93\xf8H\xa0e7\xf1\x01\xa3v\xef\x85\x95\xbd\
\xeb\n\xb7\x7f\xdc\xdd\x85\xe5\xa1If\xfa\x9b\xb9s~\xf7\xbc~\xdfs~\xdf\xf3;\xbf\
#0\xbe\xe4LM\xe7\xea\xbe_\x80\x9cb"<\xa8\xa3\x04\x05>\xf3\x99)\x9b;\xc8\x8bG\
\xd3\x87BB\x06\x9b\xd7-!;S\xc6\x9c`\x12\xd2\'\xc9B\xe9\x86]\x00\x14\x17\xce\xde\
\x7f\xb2\xdd\xb3O\xd3n_\xe0\x1e"\x8cW\xb9`&\xb7\x8a\n\xf0\xffh!\x93\xe5\x14\x12\
M&H0\x81\x0e\x84\xc2\x02b\xa2\xce\x91\x8f\x13i\xef\xb4\xa1\xde\x96\xf5\x95K\xf2\
\xf4y\xf9\xd9\xe4L\xcd\x1bRn\x85\xf5\xb7[\xdb\x87ZO\x9c3\xf5\xf5\x0f\x94\x00g\
\xbe+\x00\xc7\x96G\xf8\xa2\xbc\x88\x84\xec4D\xc1\x04\xe8`2\xc5w\xba}\x07\x06\x87\
\x04zU\xd8\xbc/\x05\xbf\x1a`\xe3\x13\xcb\x06\xb7m)\x1fLN\x16\x83\xe7\xffs\xc5\
\xbb\xf6\'\xcf;\x81U\xc0\xfb\xdf\x05@\xc1\xb2\xb9\xb8w=\x8dE\x14I\x18mx\xa4\xf4\
\x0f\xc0\xe9\xf3\x10\xba\x03\xc5\xb3\x12xr\x8f\x15\xbf\x1a\xe0\xeb\x8eCL\x92l\xdc\
\x18\x18\xc0\xe1\xaa\x00\xa8\x06\x9e\xffV\x00J\xe70\xe8t\x10\xde\xb8\x0c\x8b\xc9\
\x84 \x08\x13\x038\xf5\x05<4\x1d\x12\xcc\xa0i\x90`\x81\xae\x1b\x12M\xe7r\xd9\xf1\
\x9b\xf5,/\x9a\x8bOQ\xb8\xe0\xf97\xab7\xbd\xd0\x06\x94\x8d\x9ec\xcc\xfa\n\xa6 ,w\
\xde\xdd\xf8\xc5\xcbp\xe08\xa4\xdbA0\x81\x98\x04\xb2\x04v+\xcc\xbfO\xa5a\xddg\xac\
\xd8\xb0\x93\xd6\xb6O\xc9\x95e\xe68g\xe2>^_\xba\xda\x89~O\x00j\x10\x92\x12A\xd7\
\x8do<\xf9\xa4\x0b\x1e_\x00SeH\xb6\x8e\xdfg\xb5\x13\xd6T\xd4\xd2v\xd6M\xae,#g\xa7\
\xf1\xbb\x9a\x1aF\x83H\x18\xf9s\xff\x14\xf4\xf5\xc5\x08\x8e,\x04\xb3\xd9\x00\x10\
\xdd\x85P\x18\n\x9f\x83\xeb}\xf0\xcc\x12\x10\x13!5ub\xf7<\xe2\x82\xaf{\xe0\x0f\xfb\
\xda(+]\x88\xd3\x91\x83`\x15\x993\xdd\xc9\xf5\xee\x0fk:{\xa9\x1d\xc3\x81\x9f-C\xdf\
X\x06\xf6d\xc3\xb0\x1ea\xff\x81\xe3p\xf0}\xbe\x97d\xda\xe1z?\xe8\xde\x16\xdeln\xc3\
\xfd\xa5\x8f\xfa\xd7\x9ac\xb6c.xf)\xfa\xb4tHN\x8a\xb0S0\x8c\xff\xf5\xe4\xf77\x0e\
\x90\x96S\xc0\x91#Gx\xab\xb9\x8d\xb2%N\xb6>[Ns\xc3V\x8c\xc3=\x02\xc0\xe5\x9bP\xfa\
\x03\x83T\xd1\xd5+*\xec=:\xf1\xe4\xc9\xc9\xc9$\'\'#\x8ab\xbc\xd1\xb44222\xc8\xcc\
\xcc\xa4\xae\xae\x0e\xd7\xbc|l\x82\x82\xa0\x9b\xc8\x95e\xe6\x17\xb9b}\xcdQe\x92\x1d\
l"\x0c\rAx\x10^j\x85\xf6\x8bF\x9b\x1eac \x10\xc0n\xb7\xc7\xd5i\xda\x00CC\x90\x92\
\x92\x02@__\x1fv\xbb\x9d;\x9a\x1f\xf4A\xd0u\xbe9\xf9s\xd6V\xbeK\xcb\xdb\x0e\xa6Mv\
\xe1\xbd\xe6\x8b\x07\xf0\xe8\x83\xe8\x93eHL\x84\x81\x10<\xfa\x02\xfc\xf6\x87p\xe6\
\x12\xdc\xf8\xf2]\x00<\r\x02\xdd\xbd\xf1;\xf0\xf2O\x05\xa6\xa6\x1b\xc0\xa3b\xb7\xdb\
\xf9|_*C\xe1\xfeX\x9d\xae\x19\xa7B\xe9\xf1\x02.\x1c\xd3\x1d\xb4\x1d\xae\xa1tC\x8dn\
\x06\x90Dxj\x91\xd1\xf9\xe0qh\xddn\xe8\x8b\x0b\xc0"\x17\x90)\t\xactA\xfb\xa5QL\x7f\
`|\xd7\x8c6\x0e\xf0\xa7\n\xd8\xd1\xd8\xc8\xa6\xb5k\x00P\x14\xbf\x11\x07\x9eXL\xc7\
\xe1\x8f\x8c@\xe2\xaa\x82\xcd\xa5\xc6\x80u/\xc3\x81\xbfCj\xd6,N\xd6\xc2\xf6\xb5\xb0\
n\xf1\xf8\x06\xff\xdb\xff\xedI\xa9\x04u\xa4$\t\xe7\xac\xe9\x11\x17\x0c1\xfd\x81\\\
\xc3x\xd9\\x\xe9=x\xefcp\xd7\x1b\x03\xfe\xf8\xe4\xf0\xe0\xce\xabF\x99\x95\x95\x157\
\xa9\xefz\xbc\x11]\x03A\x1c^\xbd0\xc2E\xaa\xa6"[e|~\xc5\x00p\xae\x93[\xc10\xe9\x00\
\x1f\x9c\x1f\xee\xb8\xa3qX?\xe61\xca\xfc\x1c\xa3,,,$|k\x98\x10\xa3\xb9\x115(\xc4\
\x1f\x8e\xd8<\x88\x12\xba\xaa\xc4.#\xbd\xb9\xd2p\x81\xa2\x8e-}A\x90E\xb0\rZY\xf7j\
\x10\x80\xba\xba:f\xcf\x9e\xcd\xfe\x9d\x8fq\xcc\x03yY\xd0\xd5;|b\x9a~e\x01 73\x0c\
\x18\xe3\x01\xd4\x10\xbc\xd8l\xe8\x964\x07\xa5%.\xe3\x14\xf8\x82\xc6\x17\xd5\xf3,6\
<=\x01\xf2,6\xf2,\xd0\xa5\x05\xb0%\x0e\xaf\xa4\xa4\xa4\x04\xdb\x95\x16V\x14\xc1\x8a\
"8y6\x81\xae\xdeA\xa6M\x9bf\x18L\r\xe3\xf7\x83\xef\xba1\xc8\x07\xbc\xd3\x1e\xc6\x92\
\xe6`\x7fc\r\x12\xf0\xa1\xdb\xcd\xeb\r\r\x98\xaa\xb3]\xd0RB\xe5A8q\x16*\x0f\xc2c{\
\x03\x9c8k\x94\xfbO\x05(I\x93\t\xde0\x0fgL\x0b\x16p\xb3\xeb4rB"y\x16\x9b\x912\x01\
\xaa\xaa\xd2\xef\xbf\xcc\x9d\x80\x15\x15\xb0\x99\xcc\xb1\x0f\xe0\x8d\x03\xf5H\xa2\
\x84,K\x94\x96\x95\x1b\x1c\xe8\xc8r\xb3\xdcSN=\xe5\xe0\x81\xe5Q+Q\xdd\xa3\xd2\xe6\
\xd3\xd9\xd5w.\xce\x9fK\xb7\x9f\x8aha\xd6/2\x00\xf4\xf7\xf7#Z%D!\x81<\x8b\x8d\x81\
\xdb:\xa7\xae\xdeb\x95\xc3\x16\t\xef*R\x92\x04\x02H\xa3#\xe1\xc4"\xe1V\xaep\x13-\xae\
\xb6\x9e\xf2\x98~\xfeL\x0f`\x00\xdc\xfdb=\xd5{\x03\xb16Q\x14\xf9|\xe1C\xb8r\xbcH\xa2\
\x1c\xb9\xe7AQ\x8dhhz\xd8SzO\x08M\xfc+\xa6\x17\x17\x17\x8fi\xcf\'#\xa6WWW\xc7\xb5}\
\xd3q\x89\x9d3\x1f\x07Q@\r*\xa8!\xd5pW4!\xf9\xa7\xb3\xed\xae\xc6\xaf\xd1\x17\xf7\xbf\
{\xf7n\xaeF/\x89\x88X&\xd8\xc8\xbc\xbc<2\xe6\xe5r\xfe\xd0?\x8c\x95\x87@J\x92PCw\xc9\
\x88\xc6<\x10H\x8d\xdb\xce\xa5K\x97r\xfa\xf7\x7f\x1b\xd3o\x1b%q\xffV\xab\x95\xce\xce\
N\x94\xaek\xc3\xf7{\xe4\xe0\xebA#\x06\x14\xccvb\xce\xcf\x86*Okl\xe0\xd6r\x99=\xadJ,\
\xad\x02X\r\xecz\xe3\x13\xe6\xcf\x9fOX\x0bq\xa1\xf3#\x04Y\x88<\x14@\xd2$\xeeC\xe6\xad\
\x07\xab\xc8\xdd\xf3cf\xcc\x98\xc1\x94)S\xb8\x13\x0c\xd1\x98\xffK\xe3hf\xbb@\x94b \
\xbc\xdd^\xbe\xba\xe8\xc1\xdc\xd9\x03m\x87+q\xdc\xef\x02AB\xd7T\xaaje\xd0\x14\xf4\xa0\
\x1f\x9f\x06h*o\xbe\xb2\x8dC\x03*\xaf\xbe\xf3\xa9\x01t\x95L\xe75%\x16\xdd\xb6\xae\x04\
\xaf\x9aK\xd7\xaf\x9b\xf9\xcaU\x88\xdf\xad\x12\xf0\x86A\xd6\x91\x14\x19\x87\xe4E\x8aF\
\xc8$!\xfe:v\xcc2\x12\x84\\\xd9`\xbd\xb14\t\x90H\xd3\x00]\xc5\xe9\xa8@A\xe2\xb9Zc\xb0\
\xde\xe3%\xcd\xe1\xc0\xaf\xe8\xa8\x9a\xca\x8eg+"!\x07j\x9a\xbc\x00Tnr\xd2\xd0\x14A\xf8\
\x17X\xf3\x94\x8a\x90$\xa1\x8c\xe0\x80\x19@\x06\xa2\x91\xd9\xa7(\x110\x12\xe8 \xa1\xa2\
\x02r\x9a\x84\xac\xc72)p\xe4\xa2j~\x1cV\x01U\x84\xa3\xad-\xa0)\x90$\x19L\x0f\x1a\x19\
mM\x15\xa0\xeb\xe8\x80\x904\x92\xff#\x00<\xbd\xa9b\xec\xfbpQ)\x1dg\xda\xd8RU\x8d,\xcb\
8\xb2u\x10\x87\xb7.W\x96PB\x02\x92$CHA\xd5\x14\xa4H\xbc\x97\x92$\x10\xc7\x12Z\xd5@\xd1\
\x0c\x00\xfe^\xff\xc4\x8f\xd3\x91\xe2>\xd3\xa277\xd6\xd0q\xd6=\x16d\x91k\xdc\xfa\x91\
\xf2pi9\x15\x9b\xd7\xc6\xcfy\xb1\x9b\xd7\xff\\\xcb1\x0f\x82\xc0\xffQ\x1a_\xab\xd1\xbb\
/u\xd3\xf1A\xd3\x18\xe0\xb5\x07\xdd\x02\xc0\xff\x00+\xe4\xca\x93\x95`\xa8\xfa\x00\x00\
\x00\x00IEND\xaeB`\x82'
