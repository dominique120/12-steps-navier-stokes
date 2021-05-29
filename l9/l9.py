#!/usr/bin/env python

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.cm as cm

matplotlib.rcParams["font.family"] = "Serif"
matplotlib.rcParams["font.size"] = 10
matplotlib.rcParams["axes.labelsize"] = 10
matplotlib.rcParams["xtick.labelsize"] = 10
matplotlib.rcParams["ytick.labelsize"] = 10
matplotlib.rcParams["legend.fontsize"] = 10

fig = plt.figure(facecolor="white")
ax = fig.gca(projection='3d')
ax.grid()
ax.set_axisbelow(True)
ax.set_title("Plot of p")

x = np.array([0.0000000000000000E+00,0.3333333333333333E-01,0.6666666666666667E-01,0.1000000000000000E+00,0.1333333333333333E+00,0.1666666666666667E+00,0.2000000000000000E+00,0.2333333333333333E+00,0.2666666666666667E+00,0.3000000000000000E+00,0.3333333333333333E+00,0.3666666666666666E+00,0.4000000000000000E+00,0.4333333333333333E+00,0.4666666666666667E+00,0.5000000000000000E+00,0.5333333333333333E+00,0.5666666666666667E+00,0.6000000000000000E+00,0.6333333333333333E+00,0.6666666666666666E+00,0.7000000000000000E+00,0.7333333333333333E+00,0.7666666666666666E+00,0.8000000000000000E+00,0.8333333333333334E+00,0.8666666666666667E+00,0.9000000000000000E+00,0.9333333333333333E+00,0.9666666666666667E+00,0.1000000000000000E+01])
y = np.array([0.0000000000000000E+00,0.6666666666666667E-01,0.1333333333333333E+00,0.2000000000000000E+00,0.2666666666666667E+00,0.3333333333333333E+00,0.4000000000000000E+00,0.4666666666666667E+00,0.5333333333333333E+00,0.6000000000000000E+00,0.6666666666666666E+00,0.7333333333333333E+00,0.8000000000000000E+00,0.8666666666666667E+00,0.9333333333333333E+00,0.1000000000000000E+01,0.1066666666666667E+01,0.1133333333333333E+01,0.1200000000000000E+01,0.1266666666666667E+01,0.1333333333333333E+01,0.1400000000000000E+01,0.1466666666666667E+01,0.1533333333333333E+01,0.1600000000000000E+01,0.1666666666666667E+01,0.1733333333333333E+01,0.1800000000000000E+01,0.1866666666666667E+01,0.1933333333333333E+01,0.2000000000000000E+01])
z = np.array([np.array([0.0000000000000000E+00,0.1252395415110200E-01,0.2502260202707493E-01,0.3746548225035464E-01,0.4982088438484018E-01,0.6205525973185412E-01,0.7413261340947831E-01,0.8601386935180245E-01,0.9765619833798901E-01,0.1090122970010219E+00,0.1200296027486208E+00,0.1306494253204978E+00,0.1408059698079988E+00,0.1504252176200705E+00,0.1594236199361897E+00,0.1677065408330109E+00,0.1751663616149115E+00,0.1816801191772440E+00,0.1871064914801873E+00,0.1912818484443022E+00,0.1940149314946220E+00,0.1950794617332351E+00,0.1942035102845728E+00,0.1910535950674615E+00,0.1852097486424332E+00,0.1761241531025630E+00,0.1630475402032256E+00,0.1448863692231895E+00,0.1198947889314284E+00,0.8492367164988623E-01,0.3333333333333333E-01]),np.array([0.0000000000000000E+00,0.1252395415110200E-01,0.2502260202707493E-01,0.3746548225035464E-01,0.4982088438484018E-01,0.6205525973185412E-01,0.7413261340947831E-01,0.8601386935180245E-01,0.9765619833798901E-01,0.1090122970010219E+00,0.1200296027486208E+00,0.1306494253204978E+00,0.1408059698079988E+00,0.1504252176200705E+00,0.1594236199361897E+00,0.1677065408330109E+00,0.1751663616149115E+00,0.1816801191772440E+00,0.1871064914801873E+00,0.1912818484443022E+00,0.1940149314946220E+00,0.1950794617332351E+00,0.1942035102845728E+00,0.1910535950674615E+00,0.1852097486424332E+00,0.1761241531025630E+00,0.1630475402032256E+00,0.1448863692231895E+00,0.1198947889314284E+00,0.8492367164988623E-01,0.3333333333333333E-01]),np.array([0.0000000000000000E+00,0.1256710123079408E-01,0.2510930193572205E-01,0.3759657681992141E-01,0.4999766812275733E-01,0.6227950519943472E-01,0.7440660634817273E-01,0.8634045498886636E-01,0.9803884092295542E-01,0.1094551554693785E+00,0.1205376265890521E+00,0.1312284764024112E+00,0.1414629783621364E+00,0.1511683841600016E+00,0.1602626803145574E+00,0.1686531199263941E+00,0.1762344541750301E+00,0.1828867574366437E+00,0.1884726941190691E+00,0.1928340058470327E+00,0.1957868903931355E+00,0.1971157755310976E+00,0.1965647242190980E+00,0.1938252856932978E+00,0.1885189662066126E+00,0.1801716644867138E+00,0.1681770580528152E+00,0.1517497733136278E+00,0.1298950581025540E+00,0.1015512003413844E+00,0.6666666666666667E-01]),np.array([0.0000000000000000E+00,0.1264844212838973E-01,0.2527286039577287E-01,0.3784414176249098E-01,0.5033197173238937E-01,0.6270427063206077E-01,0.7492661398770763E-01,0.8696162999668684E-01,0.9876836599303775E-01,0.1103016142864433E+00,0.1215111856853876E+00,0.1323411162178059E+00,0.1427287887831756E+00,0.1526039463550661E+00,0.1618875664119293E+00,0.1704905568339717E+00,0.1783122206657487E+00,0.1852384197628485E+00,0.1911393441202181E+00,0.1958667636141005E+00,0.1992506027075006E+00,0.2010946430472118E+00,0.2011711482891939E+00,0.1992142950828751E+00,0.1949127008601237E+00,0.1879026728906785E+00,0.1777677469799443E+00,0.1640611302205968E+00,0.1463968394446425E+00,0.1247235296712166E+00,0.1000000000000000E+00]),np.array([0.0000000000000000E+00,0.1276704183586793E-01,0.2551139274646017E-01,0.3820531296977407E-01,0.5081990555564735E-01,0.6332457654074937E-01,0.7568648178187593E-01,0.8786995466056974E-01,0.9983590814120159E-01,0.1115412036873744E+00,0.1229379782674447E+00,0.1339729190470187E+00,0.1445864733047644E+00,0.1547119785530305E+00,0.1642746947680650E+00,0.1731907171154320E+00,0.1813657439385360E+00,0.1886936720080762E+00,0.1950549913778359E+00,0.2003149607973253E+00,0.2043215728296323E+00,0.2069033903286364E+00,0.2078675063625520E+00,0.2069982637506236E+00,0.2040582146765257E+00,0.1987946141387334E+00,0.1909584957380721E+00,0.1803505945132224E+00,0.1669199560582937E+00,0.1509501863264796E+00,0.1333333333333333E+00]),np.array([0.0000000000000000E+00,0.1292151424234884E-01,0.2582210553402513E-01,0.3867584647209535E-01,0.5145570839414858E-01,0.6413305166786426E-01,0.7667710566789430E-01,0.8905443406753082E-01,0.1012283807477111E+00,0.1131584912424690E+00,0.1247999042234869E+00,0.1361027071286874E+00,0.1470112497619475E+00,0.1574634097258421E+00,0.1673898042589693E+00,0.1767129451081335E+00,0.1853463377235363E+00,0.1931935356373770E+00,0.2001471796635214E+00,0.2060880874976570E+00,0.2108845271609504E+00,0.2143919347457011E+00,0.2164535726593003E+00,0.2169030565598349E+00,0.2155704522726662E+00,0.2122949643122135E+00,0.2069492831721621E+00,0.1994831480913261E+00,0.1899944714804747E+00,0.1788280178938605E+00,0.1666666666666667E+00]),np.array([0.0000000000000000E+00,0.1311004319297386E-01,0.2620133939713455E-01,0.3925018475335802E-01,0.5223183983718675E-01,0.6512005513222693E-01,0.7788658942001959E-01,0.9050071818921838E-01,0.1029287292368018E+00,0.1151334028197541E+00,0.1270734741059736E+00,0.1387030764105445E+00,0.1499711650175386E+00,0.1608209236509683E+00,0.1711891594813947E+00,0.1810056989424016E+00,0.1901928072505676E+00,0.1986646721476665E+00,0.2063270216324513E+00,0.2130769939715720E+00,0.2188034585795824E+00,0.2233881183137507E+00,0.2267079372363037E+00,0.2286397729721094E+00,0.2290685885459340E+00,0.2279012679920440E+00,0.2250886775666865E+00,0.2206585213428962E+00,0.2147589871140608E+00,0.2077048244412176E+00,0.2000000000000000E+00]),np.array([0.0000000000000000E+00,0.1333040934527054E-01,0.2664462369134718E-01,0.3992154104334326E-01,0.5313909727247243E-01,0.6627383062788179E-01,0.7930044233055780E-01,0.9219135175749785E-01,0.1049162451345045E+00,0.1174416173228400E+00,0.1297303073405567E+00,0.1417410299031845E+00,0.1534279078026367E+00,0.1647400138263839E+00,0.1756209368508005E+00,0.1860083957935397E+00,0.1958339389179890E+00,0.2050227870440975E+00,0.2134939112910020E+00,0.2211604845759746E+00,0.2279309189553937E+00,0.2337108078947130E+00,0.2384062438140989E+00,0.2419291802987306E+00,0.2442057338081188E+00,0.2451884837252716E+00,0.2448736978416756E+00,0.2433234839415174E+00,0.2406903144013601E+00,0.2372363571896230E+00,0.2333333333333333E+00]),np.array([0.0000000000000000E+00,0.1358002191335659E-01,0.2714674095636305E-01,0.4068199854152250E-01,0.5416675311610210E-01,0.6758068638951158E-01,0.8090180851423341E-01,0.9410606195431435E-01,0.1071669244399773E+00,0.1200550124238199E+00,0.1327376879762135E+00,0.1451786743068344E+00,0.1573376883959270E+00,0.1691701040232700E+00,0.1806266654543346E+00,0.1916532821203275E+00,0.2021909491311602E+00,0.2121758592020705E+00,0.2215398009513646E+00,0.2302109794697443E+00,0.2381154506067275E+00,0.2451794323491482E+00,0.2513328421559847E+00,0.2565144935035285E+00,0.2606794278722476E+00,0.2638087713241848E+00,0.2659221227292309E+00,0.2670915535121300E+00,0.2674545766210653E+00,0.2672210094470713E+00,0.2666666666666667E+00]),np.array([0.0000000000000000E+00,0.1385595432343563E-01,0.2770179919633413E-01,0.4152262132802984E-01,0.5530270753907780E-01,0.6902519432889191E-01,0.8267171879845993E-01,0.9622207163001036E-01,0.1096538536413496E+00,0.1229421386474682E+00,0.1360591471297094E+00,0.1489739377483647E+00,0.1616521273344105E+00,0.1740556550632893E+00,0.1861426135874290E+00,0.1978671796893149E+00,0.2091796904004390E+00,0.2200269285589912E+00,0.2303527055314780E+00,0.2400988590756858E+00,0.2492068207699919E+00,0.2576199471174231E+00,0.2652868428119250E+00,0.2721659147053521E+00,0.2782313441712478E+00,0.2834804923388609E+00,0.2879423703233186E+00,0.2916861285265351E+00,0.2948275443050988E+00,0.2975304798750170E+00,0.3000000000000000E+00]),np.array([0.0000000000000000E+00,0.1415498281373460E-01,0.2830330997147183E-01,0.4243357378254193E-01,0.5653365210309531E-01,0.7059040200698836E-01,0.8458935663503482E-01,0.9851442659152204E-01,0.1123476081171718E+00,0.1260687015276274E+00,0.1396550452056322E+00,0.1530812729614696E+00,0.1663191060772452E+00,0.1793371961424223E+00,0.1921010412604741E+00,0.2045730068136730E+00,0.2167124931931845E+00,0.2284763071354874E+00,0.2398193106600273E+00,0.2506954414701597E+00,0.2610592189902192E+00,0.2708678661899996E+00,0.2800841796415013E+00,0.2886802525314497E+00,0.2966420721784296E+00,0.3039748409861508E+00,0.3107085737393901E+00,0.3169030962566237E+00,0.3226510789491807E+00,0.3280773992047089E+00,0.3333333333333333E+00]),np.array([0.0000000000000000E+00,0.1447362708976699E-01,0.2894427046966990E-01,0.4340424562738884E-01,0.5784524019371422E-01,0.7225805184737671E-01,0.8663233060666854E-01,0.1009563272352320E+00,0.1152166502943936E+00,0.1293980355406049E+00,0.1434831330098003E+00,0.1574523194014033E+00,0.1712835464617056E+00,0.1849522401941175E+00,0.1984312711419222E+00,0.2116910229412282E+00,0.2246995949990536E+00,0.2374231854864566E+00,0.2498267124000100E+00,0.2618747419977271E+00,0.2735328028665572E+00,0.2847691655056676E+00,0.2955571538487960E+00,0.3058780150945105E+00,0.3157242926965020E+00,0.3251035093197020E+00,0.3340417653261462E+00,0.3425866123254047E+00,0.3508083377709845E+00,0.3587987301262014E+00,0.3666666666666666E+00]),np.array([0.0000000000000000E+00,0.1480819227655772E-01,0.2961724800752307E-01,0.4442338018254069E-01,0.5922226086855705E-01,0.7400880308748374E-01,0.8877694770588618E-01,0.1035194571376958E+00,0.1182277183257142E+00,0.1328915584250882E+00,0.1474990779496339E+00,0.1620365079699133E+00,0.1764881004012739E+00,0.1908360636352572E+00,0.2050605598891918E+00,0.2191397857756131E+00,0.2330501637161041E+00,0.2467666787107310E+00,0.2602634019872730E+00,0.2735142488852373E+00,0.2864940206440948E+00,0.2991797748060644E+00,0.3115525512347916E+00,0.3235994432333695E+00,0.3353159384545462E+00,0.3467083574872947E+00,0.3577960931374303E+00,0.3686132219919313E+00,0.3792089695945417E+00,0.3896465353707133E+00,0.4000000000000000E+00]),np.array([0.0000000000000000E+00,0.1515481157407923E-01,0.3031446576021331E-01,0.4547920399962846E-01,0.6064881361877805E-01,0.7582245325266049E-01,0.9099848341320897E-01,0.1061743038418151E+00,0.1213461997438938E+00,0.1365091996652132E+00,0.1516569486927211E+00,0.1667816019177618E+00,0.1818737447676032E+00,0.1969223489746018E+00,0.2119147756788916E+00,0.2268368404725621E+00,0.2416729590155538E+00,0.2564063959300680E+00,0.2710196434527318E+00,0.2854949587677886E+00,0.2998150883807439E+00,0.3139642018239055E+00,0.3279290420882710E+00,0.3417002726595165E+00,0.3552739577175481E+00,0.3686530525285104E+00,0.3818487110714393E+00,0.3948811533877481E+00,0.4077798042648708E+00,0.4205824542674723E+00,0.4333333333333333E+00]),np.array([0.0000000000000000E+00,0.1550948919975600E-01,0.3102788889530912E-01,0.4655955665820261E-01,0.6210848246685836E-01,0.7767815727138505E-01,0.9327144648199555E-01,0.1088904696909015E+00,0.1245364881461644E+00,0.1402098018048057E+00,0.1559096582307331E+00,0.1716341762056095E+00,0.1873802877260473E+00,0.2031437031054355E+00,0.2189189052091883E+00,0.2346991804259806E+00,0.2504766957561733E+00,0.2662426332247532E+00,0.2819873943834142E+00,0.2977008884042380E+00,0.3133729163226621E+00,0.3289936601364273E+00,0.3445542771738773E+00,0.3600475858122644E+00,0.3754688071922154E+00,0.3908162995666573E+00,0.4060921909273131E+00,0.4213027896280528E+00,0.4364586446312870E+00,0.4515741514749370E+00,0.4666666666666667E+00]),np.array([0.0000000000000000E+00,0.1586814337214042E-01,0.3174931064407566E-01,0.4765202009171032E-01,0.6358450869311996E-01,0.7955464358922398E-01,0.9556983806168994E-01,0.1116369729221664E+00,0.1277623241430989E+00,0.1439514974878133E+00,0.1602093708170823E+00,0.1765400446619452E+00,0.1929468015591147E+00,0.2094320745469238E+00,0.2259974251175454E+00,0.2426435308162832E+00,0.2593701825722586E+00,0.2761762917379577E+00,0.2930599067090470E+00,0.3100182388917554E+00,0.3270476976840933E+00,0.3441439340403307E+00,0.3613018920965133E+00,0.3785158682493570E+00,0.3957795770025107E+00,0.4130862228237608E+00,0.4304285771949806E+00,0.4477990599841516E+00,0.4651898242261255E+00,0.4825928433663625E+00,0.5000000000000000E+00]),np.array([0.0000000000000000E+00,0.1622664923468385E-01,0.3247043814800386E-01,0.4874404730384982E-01,0.6505996221227098E-01,0.8143042766952346E-01,0.9786740620979477E-01,0.1143825411470978E+00,0.1309871243495206E+00,0.1476920684340397E+00,0.1645078824708971E+00,0.1814446495090703E+00,0.1985120032422102E+00,0.2157190998927233E+00,0.2330745798767326E+00,0.2505865120289870E+00,0.2682623111751812E+00,0.2861086177992914E+00,0.3041311267840934E+00,0.3223343512569573E+00,0.3407213083157187E+00,0.3592931170666905E+00,0.3780485075169162E+00,0.3969832530245675E+00,0.4160895602891974E+00,0.4353554787305470E+00,0.4547644219762548E+00,0.4742949199976758E+00,0.4939207284407863E+00,0.5136113971629438E+00,0.5333333333333333E+00]),np.array([0.0000000000000000E+00,0.1658088172345669E-01,0.3318297813815987E-01,0.4982309069136170E-01,0.6651791214412106E-01,0.8328402400334912E-01,0.1001378978465403E+00,0.1170959006889884E+00,0.1341743039567118E+00,0.1513892948213007E+00,0.1687569875688970E+00,0.1862934312453997E+00,0.2040146079650182E+00,0.2219364139090605E+00,0.2400746121107601E+00,0.2584447425999302E+00,0.2770619714459323E+00,0.2959408559470894E+00,0.3150949992300286E+00,0.3345365648701150E+00,0.3542756225058513E+00,0.3743193012943288E+00,0.3946707427684581E+00,0.4153278720024546E+00,0.4362820491535964E+00,0.4575167228379409E+00,0.4790062766608170E+00,0.5007153246719093E+00,0.5225987423062640E+00,0.5446026799055824E+00,0.5666666666666667E+00]),np.array([0.0000000000000000E+00,0.1692675843767604E-01,0.3387872262008823E-01,0.5087673037502933E-01,0.6794159739504269E-01,0.8509415813256319E-01,0.1023553107249133E+00,0.1197460659378535E+00,0.1372875993238284E+00,0.1550013046883765E+00,0.1729088454594077E+00,0.1910321985509870E+00,0.2093936826769315E+00,0.2280159596574060E+00,0.2469219929356826E+00,0.2661349421853612E+00,0.2856779665522662E+00,0.3055739019747576E+00,0.3258447708008808E+00,0.3465110758815861E+00,0.3675908288025853E+00,0.3890982666826532E+00,0.4110422294917300E+00,0.4334242071890433E+00,0.4562361306152548E+00,0.4794580767403351E+00,0.5030561836117231E+00,0.5269812016703130E+00,0.5511681982784915E+00,0.5755379074290602E+00,0.6000000000000000E+00]),np.array([0.0000000000000000E+00,0.1726028257185888E-01,0.3454963473692185E-01,0.5189280294054201E-01,0.6931459806659873E-01,0.8683998018973806E-01,0.1044941479884854E+00,0.1223026328827225E+00,0.1402913970261514E+00,0.1584869329728271E+00,0.1769163610324267E+00,0.1956075178828468E+00,0.2145890267345816E+00,0.2338903350165130E+00,0.2535416999298568E+00,0.2735740950547799E+00,0.2940190023236857E+00,0.3149080431261577E+00,0.3362723905347283E+00,0.3581418928793035E+00,0.3805438297447030E+00,0.4035012196446865E+00,0.4270306119032110E+00,0.4511393351569670E+00,0.4758222562530268E+00,0.5010582412207104E+00,0.5268067111149247E+00,0.5530049314522185E+00,0.5795668974287564E+00,0.6063847435023687E+00,0.6333333333333333E+00]),np.array([0.0000000000000000E+00,0.1757758591333273E-01,0.3518793486805939E-01,0.5285953080492487E-01,0.7062100823500914E-01,0.8850128107685534E-01,0.1065296773703479E+00,0.1247360802908180E+00,0.1431510748986011E+00,0.1618060986538940E+00,0.1807335917663715E+00,0.1999671407473267E+00,0.2195416048437817E+00,0.2394932100449708E+00,0.2598595886742307E+00,0.2806797337618244E+00,0.3019938259606251E+00,0.3238428763185249E+00,0.3462681106549299E+00,0.3693100012128475E+00,0.3930068307447693E+00,0.4173926580179364E+00,0.4424945512473002E+00,0.4683289835145484E+00,0.4948973672491618E+00,0.5221808773236231E+00,0.5501350082160352E+00,0.5786847384807353E+00,0.6077216671461464E+00,0.6371048262261348E+00,0.6666666666666666E+00]),np.array([0.0000000000000000E+00,0.1787497180382714E-01,0.3578618680072095E-01,0.5376565206860644E-01,0.7184561012200068E-01,0.9005871166982188E-01,0.1084381960823711E+00,0.1270180806598361E+00,0.1458333588507217E+00,0.1649202062081100E+00,0.1843161909424932E+00,0.2040604832188794E+00,0.2241940535585262E+00,0.2447598454410550E+00,0.2658028999237818E+00,0.2873704000994111E+00,0.3095115896160646E+00,0.3322775012416385E+00,0.3557204074893340E+00,0.3798928748600502E+00,0.4048462666070030E+00,0.4306284990501808E+00,0.4572808219052770E+00,0.4848333828692289E+00,0.5132993873896993E+00,0.5426678376271420E+00,0.5728952168416501E+00,0.6038971633828217E+00,0.6355421531535840E+00,0.6676502167802071E+00,0.7000000000000000E+00]),np.array([0.0000000000000000E+00,0.1814895781897115E-01,0.3633738350743099E-01,0.5460055021979407E-01,0.7297404894837040E-01,0.9149400440513372E-01,0.1101973010250182E+00,0.1291218211531476E+00,0.1483066966855959E+00,0.1677925742705789E+00,0.1876218924931685E+00,0.2078391670416219E+00,0.2284912763690279E+00,0.2496277353619539E+00,0.2713009373508767E+00,0.2935663345091651E+00,0.3164825119710954E+00,0.3401110900430089E+00,0.3645163592918852E+00,0.3897645121434278E+00,0.4159222788737330E+00,0.4430547037545320E+00,0.4712217114386099E+00,0.5004730291105546E+00,0.5308409867932304E+00,0.5623308051413819E+00,0.5949083620201587E+00,0.6284863562473889E+00,0.6629115089637132E+00,0.6979578760011259E+00,0.7333333333333333E+00]),np.array([0.0000000000000000E+00,0.1839631774920403E-01,0.3683503170378136E-01,0.5535438247333625E-01,0.7399300689197862E-01,0.9279019536534600E-01,0.1117861622231777E+00,0.1310223323681937E+00,0.1505416469734361E+00,0.1703888901673543E+00,0.1906110374284568E+00,0.2112576245915575E+00,0.2323811336454880E+00,0.2540373874210252E+00,0.2762859391302729E+00,0.2991904334566663E+00,0.3228189018723395E+00,0.3472439335833682E+00,0.3725426312188666E+00,0.3987962115713479E+00,0.4260890386313099E+00,0.4545067689628525E+00,0.4841331380230555E+00,0.5150447168135871E+00,0.5473027426500126E+00,0.5809409636860912E+00,0.6159485684509208E+00,0.6522481979837307E+00,0.6896715918124865E+00,0.7279404324890757E+00,0.7666666666666666E+00]),np.array([0.0000000000000000E+00,0.1861412228793032E-01,0.3727323398511856E-01,0.5601820491000772E-01,0.7489037364681193E-01,0.9393184364640400E-01,0.1131857955126323E+00,0.1326968201730575E+00,0.1525112774345740E+00,0.1726776840924501E+00,0.1932471352088590E+00,0.2142737612633948E+00,0.2358152223757554E+00,0.2579332383390829E+00,0.2806941491735755E+00,0.3041694943108356E+00,0.3284365876854140E+00,0.3535790481695856E+00,0.3796872153250247E+00,0.4068583317130456E+00,0.4351962925010722E+00,0.4648106309520458E+00,0.4958141946918817E+00,0.5283186324025063E+00,0.5624263146739679E+00,0.5982166631985528E+00,0.6357242448275771E+00,0.6749060799568348E+00,0.7155981674737162E+00,0.7574695825250932E+00,0.8000000000000000E+00]),np.array([0.0000000000000000E+00,0.1879977766807468E-01,0.3764676699058524E-01,0.5658409200245656E-01,0.7565541019892300E-01,0.9490524350108313E-01,0.1143793286520753E+00,0.1341249932214540E+00,0.1541915635914090E+00,0.1746308115354339E+00,0.1954974462540782E+00,0.2168496589536251E+00,0.2387497271427749E+00,0.2612646855890037E+00,0.2844670699602626E+00,0.3084357369071353E+00,0.3332567594677380E+00,0.3590243868896760E+00,0.3858420389719952E+00,0.4138232688781159E+00,0.4430925603031826E+00,0.4737856976682324E+00,0.5060492120907570E+00,0.5400379737116816E+00,0.5759092277046064E+00,0.6138100508160665E+00,0.6538531596218495E+00,0.6960735120344397E+00,0.7403573540324252E+00,0.7863437168187923E+00,0.8333333333333334E+00]),np.array([0.0000000000000000E+00,0.1895106135738392E-01,0.3795115376217297E-01,0.5704524764256791E-01,0.7627890170530727E-01,0.9569862366170810E-01,0.1153522534311473E+00,0.1352893764175914E+00,0.1555617750565645E+00,0.1762239276604499E+00,0.1973335647013502E+00,0.2189522941268094E+00,0.2411463091762315E+00,0.2639871945217059E+00,0.2875528494268781E+00,0.3119285499244581E+00,0.3372081754152697E+00,0.3634956276456748E+00,0.3909064694580870E+00,0.4195698018935671E+00,0.4496303698232882E+00,0.4812508135839196E+00,0.5146138134737845E+00,0.5499234897139186E+00,0.5874045760229295E+00,0.6272960715153575E+00,0.6698323207683281E+00,0.7151972555845363E+00,0.7634259573699403E+00,0.8142185838090258E+00,0.8666666666666667E+00]),np.array([0.0000000000000000E+00,0.1906615384504118E-01,0.3818272829586428E-01,0.5739610448928126E-01,0.7675329488619516E-01,0.9630232750401681E-01,0.1160926552068153E+00,0.1361755988687037E+00,0.1566048341943444E+00,0.1774369333532272E+00,0.1987319744899755E+00,0.2205542347931606E+00,0.2429729863950658E+00,0.2660634192903606E+00,0.2899077221956722E+00,0.3145963614992814E+00,0.3402296110725160E+00,0.3669194028790754E+00,0.3947915913243201E+00,0.4239887541410788E+00,0.4546736885857633E+00,0.4870337971440970E+00,0.5212865673512158E+00,0.5576862605638661E+00,0.5965315176484061E+00,0.6381722555894713E+00,0.6830102849848051E+00,0.7314770322497969E+00,0.7839425730335574E+00,0.8404419806339498E+00,0.9000000000000000E+00]),np.array([0.0000000000000000E+00,0.1914366554241929E-01,0.3833869025390497E-01,0.5763240839067851E-01,0.7707281520568654E-01,0.9670896745362971E-01,0.1165914108272135E+00,0.1367726440355562E+00,0.1573076304822127E+00,0.1782543702876911E+00,0.1996745479203088E+00,0.2216342742294851E+00,0.2442049466281439E+00,0.2674642581542148E+00,0.2914973960657386E+00,0.3163984848708538E+00,0.3422723493902755E+00,0.3692367039275773E+00,0.3974249191898594E+00,0.4269895878532693E+00,0.4581072167000008E+00,0.4909845412755611E+00,0.5258672255975323E+00,0.5630521311011457E+00,0.6029049796927700E+00,0.6458860643003691E+00,0.6925870191188511E+00,0.7437778148587384E+00,0.8004372584861678E+00,0.8636107518411177E+00,0.9333333333333333E+00]),np.array([0.0000000000000000E+00,0.1918265788766945E-01,0.3841714794063884E-01,0.5775128483259779E-01,0.7723355937466826E-01,0.9691354734504307E-01,0.1168423462353133E+00,0.1370730500984653E+00,0.1576612740801526E+00,0.1786657422270891E+00,0.2001489550199120E+00,0.2221779571080677E+00,0.2448252314320865E+00,0.2681697540603320E+00,0.2922982558023854E+00,0.3173067538013568E+00,0.3433024417470752E+00,0.3704060658358027E+00,0.3987549731369424E+00,0.4285071135819049E+00,0.4598464316687260E+00,0.4929903470255675E+00,0.5282004892252539E+00,0.5657987213590431E+00,0.6061922063237496E+00,0.6499149183130535E+00,0.6977013996508220E+00,0.7506297487003091E+00,0.8104298306141640E+00,0.8802344210082420E+00,0.9666666666666667E+00]),np.array([0.0000000000000000E+00,0.1918265788766945E-01,0.3841714794063884E-01,0.5775128483259779E-01,0.7723355937466826E-01,0.9691354734504307E-01,0.1168423462353133E+00,0.1370730500984653E+00,0.1576612740801526E+00,0.1786657422270891E+00,0.2001489550199120E+00,0.2221779571080677E+00,0.2448252314320865E+00,0.2681697540603320E+00,0.2922982558023854E+00,0.3173067538013568E+00,0.3433024417470752E+00,0.3704060658358027E+00,0.3987549731369424E+00,0.4285071135819049E+00,0.4598464316687260E+00,0.4929903470255675E+00,0.5282004892252539E+00,0.5657987213590431E+00,0.6061922063237496E+00,0.6499149183130535E+00,0.6977013996508220E+00,0.7506297487003091E+00,0.8104298306141640E+00,0.8802344210082420E+00,0.9666666666666667E+00])])

X, Y = np.meshgrid(x, y)
Z = np.transpose(z)
CS = ax.plot_surface(X,Y,Z,label="value of p",antialiased=False,cmap=cm.cividis)

plt.savefig("lp.png", dpi=200)

