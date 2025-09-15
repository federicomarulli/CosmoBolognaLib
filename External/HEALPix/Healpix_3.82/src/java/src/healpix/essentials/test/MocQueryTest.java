/*
 *  This file is part of Healpix Java.
 *
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This code is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this code; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.sourceforge.net
 */
package healpix.essentials.test;

import healpix.essentials.*;
import junit.framework.TestCase;
import java.io.*;
import java.util.Arrays;
import java.util.ArrayList;

/** Tests for the MocQuery class

    @copyright 2015 Max-Planck-Society
    @author Martin Reinecke */
public class MocQueryTest extends TestCase
  {
  private ArrayList<Vec3> rawToPoly(double[] inp) throws Exception
    {
    HealpixUtils.check((inp.length>=6)&&((inp.length&1)==0),
      "bad input array");
    ArrayList<Vec3> res = new ArrayList<Vec3>();
    for (int i=0; i<inp.length; i+=2)
      res.add(new Vec3(new Pointing(Math.PI/2 - Math.toRadians(inp[i+1]),
        Math.toRadians(inp[i]))));
    return res;
    }
  private void checkPoly(ArrayList<Vec3> poly, int o1, int o2, String out)
    throws Exception
    {
    Moc moc = MocQuery.queryGeneralPolygonInclusive (poly,o1,o1);
    if (out!="")
      MocFitsIO.mocToFits(moc,out+"_inc_"+o1+".fits");
    for (int o=o1+1; o<=o2; ++o)
      {
      Moc moc2 = MocQuery.queryGeneralPolygonInclusive (poly,o1,o2);
      assertTrue("inconsistency",moc.contains(moc2));
      if (out!="")
        MocFitsIO.mocToFits(moc2,out+"_inc_"+o+".fits");
      moc=moc2;
      }
    Moc moc2 = MocQuery.queryGeneralPolygon(poly,o1);
    if (out!="")
      MocFitsIO.mocToFits(moc2,out+".fits");
    assertTrue("inconsistency",moc.contains(moc2));
    }
  public void test1() throws Exception
    {
    double[] polyraw=new double[]{
      83.6479647404001,  22.06796566083695,
      83.57919759304697, 22.014902754311912,
      83.54565527259793, 22.05666932013043,
      83.54535325267885, 22.05666728060237,
      83.54535545323692, 22.056387364445676,
      83.54535985432697, 22.05582753212609,
      83.54536865640274, 22.05470786746255,
      83.54537525786826, 22.05386811894398,
      83.53952949498002, 21.952773396955234,
      83.63157603991105, 21.995362438376056,
      83.64144959973622, 21.965751089665403,
      83.69912966373616, 22.00612703969729 };

    checkPoly(rawToPoly(polyraw),13,15,"");
    }
  public void test2() throws Exception
    {
    double[] polyraw=new double[]{
      83.69572488226052, 22.01842407483951,
      83.66469321144844, 22.05155150183524,
      83.66439120028902, 22.051549680680765,
      83.62952475535761, 22.02838142338798,
      83.58771024980035, 22.04715033027345,
      83.5874082491367,  22.04714836793211,
      83.58741036646086, 22.0468684509572,
      83.58741671838312, 22.046028700022926,
      83.57147626994606, 21.99805614401068,
      83.59369254263831, 21.934377044862455,
      83.59369043939574, 21.93465696148893,
      83.59368833614487, 21.934936878119167,
      83.59338236664615, 21.935494760357145,
      83.59338026277467, 21.93577467699384,
      83.64687451268995, 21.966904215256918,
      83.60550939897153, 22.00919444338953,
      83.66653272555045, 22.004534367269507,
      83.69542484507058, 22.018142392123163 };

    checkPoly(rawToPoly(polyraw),13,15,"");
    }

  public void test3() throws Exception
    {
    ArrayList<Vec3> vert=new ArrayList<Vec3>();
    vert.add(new Vec3(new Pointing(1.1523428941043317,1.000795791041251)));
    vert.add(new Vec3(new Pointing(1.1523106929911657,1.0008281306683644)));
    vert.add(new Vec3(new Pointing(1.1521818924889993,1.0009574985526972)));
    vert.add(new Vec3(new Pointing(1.1519887035905163,1.001151578519077)));
    vert.add(new Vec3(new Pointing(1.151859918897493,1.0012809839336458)));
    vert.add(new Vec3(new Pointing(1.1514894198468866,1.001263637028638)));
    vert.add(new Vec3(new Pointing(1.1513659202048774,1.0012578536443903)));
    vert.add(new Vec3(new Pointing(1.1511859573237597,1.0011139425394953)));
    vert.add(new Vec3(new Pointing(1.1510060020113047,1.0009700085527964)));
    vert.add(new Vec3(new Pointing(1.1508556031353578,1.0008613172356937)));
    vert.add(new Vec3(new Pointing(1.150673007120061,1.0007849744501682)));
    vert.add(new Vec3(new Pointing(1.1503938085112773,1.0008057187152946)));
    vert.add(new Vec3(new Pointing(1.150205908484378,1.0008646491036832)));
    vert.add(new Vec3(new Pointing(1.149953611594119,1.000988341677547)));
    vert.add(new Vec3(new Pointing(1.1497335179693469,1.0010796806589697)));
    vert.add(new Vec3(new Pointing(1.1494247698862567,1.0010651879076278)));
    vert.add(new Vec3(new Pointing(1.1492448165599576,1.0009211385295873)));
    vert.add(new Vec3(new Pointing(1.149000469199188,1.0008418420185987)));
    vert.add(new Vec3(new Pointing(1.1485894925444542,1.000279935669881)));
    vert.add(new Vec3(new Pointing(1.1486539077047668,1.0002151640462313)));
    vert.add(new Vec3(new Pointing(1.1490324086471566,1.0000296562389843)));
    vert.add(new Vec3(new Pointing(1.1492203271505665,0.9999707639528252)));
    vert.add(new Vec3(new Pointing(1.1495639559481354,0.9998853760122856)));
    vert.add(new Vec3(new Pointing(1.149848517849616,0.9997294255486839)));
    vert.add(new Vec3(new Pointing(1.1499773768285775,0.9995999921162992)));
    vert.add(new Vec3(new Pointing(1.1501062421175312,0.9994705737522844)));
    vert.add(new Vec3(new Pointing(1.1501733657022064,0.9993382245815566)));
    vert.add(new Vec3(new Pointing(1.150157309455921,0.9989646985587307)));
    vert.add(new Vec3(new Pointing(1.1502783613271745,0.9982588937208289)));
    vert.add(new Vec3(new Pointing(1.150367358658019,0.9975854875791592)));
    vert.add(new Vec3(new Pointing(1.1503145561122339,0.9966002329829549)));
    vert.add(new Vec3(new Pointing(1.1501837017642849,0.9960177261452899)));
    vert.add(new Vec3(new Pointing(1.1500880289072624,0.9953352142578803)));
    vert.add(new Vec3(new Pointing(1.1498339618638964,0.9947464292684306)));
    vert.add(new Vec3(new Pointing(1.1496599663238942,0.9944665891792631)));
    vert.add(new Vec3(new Pointing(1.1495211017074738,0.9940867625529576)));
    vert.add(new Vec3(new Pointing(1.1494032633666398,0.993945254823653)));
    vert.add(new Vec3(new Pointing(1.1492854324600088,0.9938037321908275)));
    vert.add(new Vec3(new Pointing(1.1491647775622915,0.9937298631728655)));
    vert.add(new Vec3(new Pointing(1.1489093237140704,0.9939204821489465)));
    vert.add(new Vec3(new Pointing(1.1487240963858985,0.9939111896936369)));
    vert.add(new Vec3(new Pointing(1.148418211419246,0.9938280074878991)));
    vert.add(new Vec3(new Pointing(1.148203527561705,0.9937833060402974)));
    vert.add(new Vec3(new Pointing(1.1480562505373242,0.9936062948138628)));
    vert.add(new Vec3(new Pointing(1.1479678899060295,0.9935000768652547)));
    vert.add(new Vec3(new Pointing(1.1479089851456499,0.9934292602260698)));
    vert.add(new Vec3(new Pointing(1.1478263164713338,0.9931876107503619)));
    vert.add(new Vec3(new Pointing(1.147658188833066,0.9927719711268423)));
    vert.add(new Vec3(new Pointing(1.147605019627341,0.9925656960713723)));
    vert.add(new Vec3(new Pointing(1.147513838399067,0.9925271325090353)));
    vert.add(new Vec3(new Pointing(1.1469790356795493,0.9927375232998785)));
    vert.add(new Vec3(new Pointing(1.1467586597155184,0.9928281691313717)));
    vert.add(new Vec3(new Pointing(1.1466646212864107,0.9928573452285496)));
    vert.add(new Vec3(new Pointing(1.146570583168876,0.9928865238140747)));
    vert.add(new Vec3(new Pointing(1.1464822508174166,0.9927802065018877)));
    vert.add(new Vec3(new Pointing(1.1464233649158908,0.9927093236024059)));
    vert.add(new Vec3(new Pointing(1.1463350395636842,0.9926029922143758)));
    vert.add(new Vec3(new Pointing(1.1464082297063476,0.9923351250888106)));
    vert.add(new Vec3(new Pointing(1.1465403155639036,0.9921381737414078)));
    vert.add(new Vec3(new Pointing(1.146801650920068,0.9918120793817878)));
    vert.add(new Vec3(new Pointing(1.146992632145478,0.9916860897831967)));
    vert.add(new Vec3(new Pointing(1.147336526500701,0.9916018945040467)));
    vert.add(new Vec3(new Pointing(1.1475540582351127,0.991579108720834)));
    vert.add(new Vec3(new Pointing(1.148297344412475,0.9908370227283142)));
    vert.add(new Vec3(new Pointing(1.14831771135701,0.9903631889078771)));
    vert.add(new Vec3(new Pointing(1.1483410904887272,0.9898216751345671)));
    vert.add(new Vec3(new Pointing(1.14817642568152,0.9893381849506938)));
    vert.add(new Vec3(new Pointing(1.1478619397258643,0.9887449725688633)));
    vert.add(new Vec3(new Pointing(1.1476621504872502,0.9883612311158905)));
    vert.add(new Vec3(new Pointing(1.1473389454791498,0.987970892009681)));
    vert.add(new Vec3(new Pointing(1.1470745474697959,0.9876514390156884)));
    vert.add(new Vec3(new Pointing(1.1467220762850698,0.9872253830842613)));
    vert.add(new Vec3(new Pointing(1.1465134990729606,0.9870445283107344)));
    vert.add(new Vec3(new Pointing(1.1463576636239432,0.9870701715700294)));
    vert.add(new Vec3(new Pointing(1.1461311065061512,0.987295784658764)));
    vert.add(new Vec3(new Pointing(1.145969291549504,0.9874569649988388)));
    vert.add(new Vec3(new Pointing(1.1456780489082319,0.9877471487418389)));
    vert.add(new Vec3(new Pointing(1.1453544826095992,0.9880696643744283)));
    vert.add(new Vec3(new Pointing(1.1451309779989192,0.9882276968098427)));
    vert.add(new Vec3(new Pointing(1.1449223325837785,0.9880468062013319)));
    vert.add(new Vec3(new Pointing(1.144684321314146,0.9878303468533963)));
    vert.add(new Vec3(new Pointing(1.1444786880119573,0.9875815784705101)));
    vert.add(new Vec3(new Pointing(1.1444055321874065,0.9871359210384268)));
    vert.add(new Vec3(new Pointing(1.1443527978610553,0.9869292081606309)));
    vert.add(new Vec3(new Pointing(1.1443588066568586,0.9867935954286845)));
    vert.add(new Vec3(new Pointing(1.1443294428641322,0.9867580408743656)));
    vert.add(new Vec3(new Pointing(1.1441178863108807,0.9866447588708069)));
    vert.add(new Vec3(new Pointing(1.143627025593574,0.9865504721602422)));
    vert.add(new Vec3(new Pointing(1.1430714198687517,0.9865206748900628)));
    vert.add(new Vec3(new Pointing(1.142486453971728,0.9864552736851653)));
    vert.add(new Vec3(new Pointing(1.1421130353865951,0.9865032583536448)));
    vert.add(new Vec3(new Pointing(1.1418013514800418,0.9865545719393494)));
    vert.add(new Vec3(new Pointing(1.1415280685012117,0.9864378231638371)));
    vert.add(new Vec3(new Pointing(1.1413195436592867,0.9862564662025781)));
    vert.add(new Vec3(new Pointing(1.1411140566169835,0.9860071695710151)));
    vert.add(new Vec3(new Pointing(1.1410583793726072,0.9858680266790968)));
    vert.add(new Vec3(new Pointing(1.1411585594318707,0.9857032675075046)));
    vert.add(new Vec3(new Pointing(1.1415471911174158,0.9853159421723982)));
    vert.add(new Vec3(new Pointing(1.1417091373202277,0.9851545968685336)));
    vert.add(new Vec3(new Pointing(1.141935878135893,0.9849287531836823)));
    vert.add(new Vec3(new Pointing(1.1421919741278386,0.9847385729150904)));
    vert.add(new Vec3(new Pointing(1.1427212823451636,0.9846653708690231)));
    vert.add(new Vec3(new Pointing(1.1429388730086534,0.9846432208437782)));
    vert.add(new Vec3(new Pointing(1.1431888603402285,0.9845888377168718)));
    vert.add(new Vec3(new Pointing(1.143221257341095,0.9845566013231962)));
    vert.add(new Vec3(new Pointing(1.1431393920272823,0.9843141087835374)));
    vert.add(new Vec3(new Pointing(1.142940242012983,0.9839291446460309)));
    vert.add(new Vec3(new Pointing(1.1427704685698539,0.9835797337470248)));
    vert.add(new Vec3(new Pointing(1.1426917873195213,0.9832693018292363)));
    vert.add(new Vec3(new Pointing(1.142613142293535,0.9829588475331977)));
    vert.add(new Vec3(new Pointing(1.1426935245487904,0.9825551447276581)));
    vert.add(new Vec3(new Pointing(1.1428587549013667,0.9823262384230973)));
    vert.add(new Vec3(new Pointing(1.1432802997807938,0.9819076068355195)));
    vert.add(new Vec3(new Pointing(1.1435428818201463,0.981582241584326)));
    vert.add(new Vec3(new Pointing(1.1438023569877103,0.9813247693528765)));
    vert.add(new Vec3(new Pointing(1.1440879942426359,0.981170799141964)));
    vert.add(new Vec3(new Pointing(1.1443087617014966,0.9810812052364875)));
    vert.add(new Vec3(new Pointing(1.1444939405572025,0.9810915890969772)));
    vert.add(new Vec3(new Pointing(1.1446435321631656,0.9812019262076423)));
    vert.add(new Vec3(new Pointing(1.1449657419610033,0.9815937325674373)));
    vert.add(new Vec3(new Pointing(1.145170814762998,0.9818430040096161)));
    vert.add(new Vec3(new Pointing(1.145405212420505,0.9821278286770655)));
    vert.add(new Vec3(new Pointing(1.1455810307686392,0.982341407316913)));
    vert.add(new Vec3(new Pointing(1.1456396407122253,0.9824125926080799)));
    vert.add(new Vec3(new Pointing(1.1457568663384885,0.9825549518114552)));
    vert.add(new Vec3(new Pointing(1.1458447905784235,0.982661711260155)));
    vert.add(new Vec3(new Pointing(1.146056176224855,0.9827752848212115)));
    vert.add(new Vec3(new Pointing(1.1464203363305736,0.982931223860804)));
    vert.add(new Vec3(new Pointing(1.1465700046565381,0.983041338988599)));
    vert.add(new Vec3(new Pointing(1.146931084242241,0.9832649203344934)));
    vert.add(new Vec3(new Pointing(1.147304554977569,0.983217580321881)));
    vert.add(new Vec3(new Pointing(1.14755456833011,0.983163460791517)));
    vert.add(new Vec3(new Pointing(1.1476518102952684,0.9830670054161794)));
    vert.add(new Vec3(new Pointing(1.147911138996274,0.9828098319777566)));
    vert.add(new Vec3(new Pointing(1.1480146153175843,0.9825780284118409)));
    vert.add(new Vec3(new Pointing(1.1478512630860116,0.982093928824457)));
    vert.add(new Vec3(new Pointing(1.1477465895705665,0.981680894029768)));
    vert.add(new Vec3(new Pointing(1.1477037064890643,0.9812712698299393)));
    vert.add(new Vec3(new Pointing(1.1476901713851657,0.980897208503522)));
    vert.add(new Vec3(new Pointing(1.147738413955949,0.9805266127850084)));
    vert.add(new Vec3(new Pointing(1.1479718818582025,0.9801664735990022)));
    vert.add(new Vec3(new Pointing(1.148205397839398,0.9798064091969002)));
    vert.add(new Vec3(new Pointing(1.1483708576902831,0.9795782695316553)));
    vert.add(new Vec3(new Pointing(1.1487928200059858,0.9791611823574279)));
    vert.add(new Vec3(new Pointing(1.1493382912635042,0.9787512843635041)));
    vert.add(new Vec3(new Pointing(1.1502609222180684,0.9788716216588957)));
    vert.add(new Vec3(new Pointing(1.1504396821178993,0.9790173766494347)));
    vert.add(new Vec3(new Pointing(1.1507908399518807,0.9794440181022032)));
    vert.add(new Vec3(new Pointing(1.1509893471305994,0.9798280210118206)));
    vert.add(new Vec3(new Pointing(1.1511815695494167,0.9803471360053592)));
    vert.add(new Vec3(new Pointing(1.1513477656945341,0.9807630544316528)));
    vert.add(new Vec3(new Pointing(1.1516604798024512,0.98135646536512)));
    vert.add(new Vec3(new Pointing(1.1517033801181038,0.981765361770002)));
    vert.add(new Vec3(new Pointing(1.1517170406506425,0.9821387451013951)));
    vert.add(new Vec3(new Pointing(1.1515842233393097,0.9823346503008543)));
    vert.add(new Vec3(new Pointing(1.1513896921866162,0.9825271705205277)));
    vert.add(new Vec3(new Pointing(1.1511627565046656,0.9827518193104232)));
    vert.add(new Vec3(new Pointing(1.1508710094941494,0.9830407197788995)));
    vert.add(new Vec3(new Pointing(1.1506441165054242,0.9832654717657893)));
    vert.add(new Vec3(new Pointing(1.150482061529942,0.9834260365777083)));
    vert.add(new Vec3(new Pointing(1.1503524244193444,0.9835545050450514)));
    vert.add(new Vec3(new Pointing(1.1502258785647448,0.983615364692989)));
    vert.add(new Vec3(new Pointing(1.150040689827067,0.9836052230787131)));
    vert.add(new Vec3(new Pointing(1.1497289561324662,0.9836559547866632)));
    vert.add(new Vec3(new Pointing(1.1496348192975478,0.9836847052609023)));
    vert.add(new Vec3(new Pointing(1.1494727890458836,0.9838453720714817)));
    vert.add(new Vec3(new Pointing(1.1494049018076335,0.9839772939956009)));
    vert.add(new Vec3(new Pointing(1.149324740233309,0.9843798294222985)));
    vert.add(new Vec3(new Pointing(1.1492739766976383,0.9848178951560865)));
    vert.add(new Vec3(new Pointing(1.1492880620951265,0.9851916630592842)));
    vert.add(new Vec3(new Pointing(1.1494023951081926,0.9854012904878543)));
    vert.add(new Vec3(new Pointing(1.149519781095103,0.9855432489469033)));
    vert.add(new Vec3(new Pointing(1.1497839272856774,0.9858626006164312)));
    vert.add(new Vec3(new Pointing(1.1501482692615095,0.986017767907073)));
    vert.add(new Vec3(new Pointing(1.1506158045877952,0.9859411952965156)));
    vert.add(new Vec3(new Pointing(1.1510863719564317,0.9857970535978254)));
    vert.add(new Vec3(new Pointing(1.1514010989298469,0.9856784621213586)));
    vert.add(new Vec3(new Pointing(1.1517158310959092,0.9855599035665571)));
    vert.add(new Vec3(new Pointing(1.1520923005998371,0.9854447044733927)));
    vert.add(new Vec3(new Pointing(1.152556811551436,0.9854358697505938)));
    vert.add(new Vec3(new Pointing(1.152936329868326,0.9852531893604244)));
    vert.add(new Vec3(new Pointing(1.1534980080389958,0.9851480997557343)));
    vert.add(new Vec3(new Pointing(1.1542448858126413,0.9850530625072047)));
    vert.add(new Vec3(new Pointing(1.1548713566503563,0.9848839302619462)));
    vert.add(new Vec3(new Pointing(1.155238694559804,0.9849714230186063)));
    vert.add(new Vec3(new Pointing(1.155511913488695,0.9850876232716619)));
    vert.add(new Vec3(new Pointing(1.1556879733525678,0.9852999969140172)));
    vert.add(new Vec3(new Pointing(1.1559784112498703,0.9857213321874994)));
    vert.add(new Vec3(new Pointing(1.1559865844022958,0.9862288863472906)));
    vert.add(new Vec3(new Pointing(1.1560385791492556,0.9871444375227638)));
    vert.add(new Vec3(new Pointing(1.1560116805212006,0.9877515244318255)));
    vert.add(new Vec3(new Pointing(1.1558379990039347,0.9881818725706057)));
    vert.add(new Vec3(new Pointing(1.155543873841697,0.9885383490226638)));
    vert.add(new Vec3(new Pointing(1.1552557042400238,0.9887599602232404)));
    vert.add(new Vec3(new Pointing(1.1549352096290768,0.9890137597369755)));
    vert.add(new Vec3(new Pointing(1.1546147391191914,0.9892676311236104)));
    vert.add(new Vec3(new Pointing(1.1542354927708685,0.9894508532485825)));
    vert.add(new Vec3(new Pointing(1.1538621301597058,0.9894990936671686)));
    vert.add(new Vec3(new Pointing(1.1536181187451178,0.9894187147674598)));
    vert.add(new Vec3(new Pointing(1.1533447121031735,0.9893029412643478)));
    vert.add(new Vec3(new Pointing(1.1529537305101472,0.9890455986187472)));
    vert.add(new Vec3(new Pointing(1.1524980830083356,0.9888524937214215)));
    vert.add(new Vec3(new Pointing(1.1523481741773312,0.9887430639190622)));
    vert.add(new Vec3(new Pointing(1.1520629738824093,0.988897439436812)));
    vert.add(new Vec3(new Pointing(1.1518395192779805,0.9890550816485132)));
    vert.add(new Vec3(new Pointing(1.1515513932259047,0.9892771221299574)));
    vert.add(new Vec3(new Pointing(1.1513250228779512,0.9895024356245301)));
    vert.add(new Vec3(new Pointing(1.1510957407699558,0.9897953995250245)));
    vert.add(new Vec3(new Pointing(1.1508988210950766,0.9900562178915127)));
    vert.add(new Vec3(new Pointing(1.150763665085474,0.9903202764372446)));
    vert.add(new Vec3(new Pointing(1.1507843367123325,0.990558526439852)));
    vert.add(new Vec3(new Pointing(1.1507992237525468,0.9909320025905168)));
    vert.add(new Vec3(new Pointing(1.1508522671724446,0.9911380213792874)));
    vert.add(new Vec3(new Pointing(1.1509376417668262,0.9913118073666528)));
    vert.add(new Vec3(new Pointing(1.1510230276226834,0.9914855800960104)));
    vert.add(new Vec3(new Pointing(1.1511701646336845,0.9916624962872871)));
    vert.add(new Vec3(new Pointing(1.1513467444499845,0.9918747650352789)));
    vert.add(new Vec3(new Pointing(1.151493907126133,0.9920516301034698)));
    vert.add(new Vec3(new Pointing(1.15167051772407,0.9922638375354705)));
    vert.add(new Vec3(new Pointing(1.1517825427272692,0.9925404659613722)));
    vert.add(new Vec3(new Pointing(1.1518945962236837,0.9928170667717601)));
    vert.add(new Vec3(new Pointing(1.151944936626454,0.9930905238387734)));
    vert.add(new Vec3(new Pointing(1.1519364038488704,0.9932932756669433)));
    vert.add(new Vec3(new Pointing(1.1519545044629247,0.9935989574717053)));
    vert.add(new Vec3(new Pointing(1.1518899290891478,0.9936634411834805)));
    vert.add(new Vec3(new Pointing(1.1518576419864373,0.9936956844356148)));
    vert.add(new Vec3(new Pointing(1.1515670755926695,0.9939859156123069)));
    vert.add(new Vec3(new Pointing(1.1514056634584344,0.9941471877691483)));
    vert.add(new Vec3(new Pointing(1.151209167317511,0.9944083526588667)));
    vert.add(new Vec3(new Pointing(1.1510098905412778,0.9947371775119466)));
    vert.add(new Vec3(new Pointing(1.1509074670044037,0.9949692293531077)));
    vert.add(new Vec3(new Pointing(1.1507939027845748,0.9954717874109883)));
    vert.add(new Vec3(new Pointing(1.1507154681212395,0.9958744790360815)));
    vert.add(new Vec3(new Pointing(1.150728328381526,0.9963155562525036)));
    vert.add(new Vec3(new Pointing(1.1507145219246022,0.9966536818156724)));
    vert.add(new Vec3(new Pointing(1.150698010346289,0.9970594382294491)));
    vert.add(new Vec3(new Pointing(1.1507515235186312,0.9972653214392541)));
    vert.add(new Vec3(new Pointing(1.15074330622707,0.9974681979653209)));
    vert.add(new Vec3(new Pointing(1.150794120124925,0.9977416918927585)));
    vert.add(new Vec3(new Pointing(1.1508449617882088,0.9980151735411477)));
    vert.add(new Vec3(new Pointing(1.150810006084698,0.9981151193696918)));
    vert.add(new Vec3(new Pointing(1.1508045689258863,0.9982503689627922)));
    vert.add(new Vec3(new Pointing(1.1507669084789147,0.9984179465727627)));
    vert.add(new Vec3(new Pointing(1.1507882982747644,0.9986561219049502)));
    vert.add(new Vec3(new Pointing(1.1508097091375338,0.9988942928018621)));
    vert.add(new Vec3(new Pointing(1.1508955843794588,0.9990677867461779)));
    vert.add(new Vec3(new Pointing(1.1509814707456794,0.9992412674747779)));
    vert.add(new Vec3(new Pointing(1.1510968986908656,0.999450015299713)));
    vert.add(new Vec3(new Pointing(1.1511559609826625,0.9995205731460943)));
    vert.add(new Vec3(new Pointing(1.151370737557223,0.9995646668380104)));
    vert.add(new Vec3(new Pointing(1.1514620182432866,0.9996028778795772)));
    vert.add(new Vec3(new Pointing(1.1515855148424001,0.9996087523184793)));
    vert.add(new Vec3(new Pointing(1.1516794788958802,0.9995793575868965)));
    vert.add(new Vec3(new Pointing(1.15174122717978,0.9995822953104395)));
    vert.add(new Vec3(new Pointing(1.1518647237644064,0.9995881703312073)));
    vert.add(new Vec3(new Pointing(1.151894256435611,0.9996234359315822)));
    vert.add(new Vec3(new Pointing(1.1519855379350588,0.9996616359063851)));
    vert.add(new Vec3(new Pointing(1.152015071556617,0.9996968986628608)));
    checkPoly(vert,13,17,"");
    }
  }
