import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;




public class TestReadNewick {
	
	public int numberOfTaxa;

	public static void main(String[] args){

		TestReadNewick trn = new TestReadNewick();
		int[][] tree = trn.makeNewickTree("");
		trn.numberOfTaxa = 8;
		
		ArrayList<Integer> members = trn.getMembers(tree, 0);
		
		for(int i = 0; i < tree.length; i++){
			System.out.println("" + tree[i][0] + " " + tree[i][1] + " " + tree[i][2]);
		}
		
		System.out.println("---");
		
		for(int i = 0; i < members.size(); i++){
			
		}
	}


	public int[][] makeNewickTree(String input){
//		input = "((((((27:1002.805899[&&NHX:age=0.000000],(86:0.000000[&&NHX:age=0.000000],13:0.000000[&&NHX:age=0.000000])172:1002.805899[&&NHX:age=0.000000])107:0.000000[&&NHX:age=1002.805899],1:1002.805899[&&NHX:age=0.000000])101:0.000000[&&NHX:age=1002.805899],((63:0.000000[&&NHX:age=0.000000],32:0.000000[&&NHX:age=0.000000])126:0.000000[&&NHX:age=0.000000],(((94:0.000000[&&NHX:age=0.000000],((90:0.000000[&&NHX:age=0.000000],67:0.000000[&&NHX:age=0.000000])180:0.000000[&&NHX:age=0.000000],(51:0.000000[&&NHX:age=0.000000],40:0.000000[&&NHX:age=0.000000])102:0.000000[&&NHX:age=0.000000])134:0.000000[&&NHX:age=0.000000])188:0.000000[&&NHX:age=0.000000],(43:0.000000[&&NHX:age=0.000000],((((99:0.000000[&&NHX:age=0.000000],70:0.000000[&&NHX:age=0.000000])198:0.000000[&&NHX:age=0.000000],37:0.000000[&&NHX:age=0.000000])140:0.000000[&&NHX:age=0.000000],(85:0.000000[&&NHX:age=0.000000],55:0.000000[&&NHX:age=0.000000])110:0.000000[&&NHX:age=0.000000])170:0.000000[&&NHX:age=0.000000],36:0.000000[&&NHX:age=0.000000])147:0.000000[&&NHX:age=0.000000])171:0.000000[&&NHX:age=0.000000])159:0.000000[&&NHX:age=0.000000],((60:0.000000[&&NHX:age=0.000000],(78:0.000000[&&NHX:age=0.000000],(96:0.000000[&&NHX:age=0.000000],(30:0.000000[&&NHX:age=0.000000],21:0.000000[&&NHX:age=0.000000])165:0.000000[&&NHX:age=0.000000])192:0.000000[&&NHX:age=0.000000])156:0.000000[&&NHX:age=0.000000])120:0.000000[&&NHX:age=0.000000],7:0.000000[&&NHX:age=0.000000])169:0.000000[&&NHX:age=0.000000])143:0.000000[&&NHX:age=0.000000])127:1002.805899[&&NHX:age=0.000000])158:4361.040322[&&NHX:age=1002.805899],11:5363.846221[&&NHX:age=0.000000])119:2687.856147[&&NHX:age=5363.846221],((52:1002.805899[&&NHX:age=0.000000],(54:639.178343[&&NHX:age=0.000000],92:639.178343[&&NHX:age=0.000000])108:363.627557[&&NHX:age=639.178343])133:2559.449441[&&NHX:age=1002.805899],(22:1002.805899[&&NHX:age=0.000000],(19:1002.805899[&&NHX:age=0.000000],33:1002.805899[&&NHX:age=0.000000])149:0.000000[&&NHX:age=1002.805899])173:2559.449441[&&NHX:age=1002.805899])160:4489.447028[&&NHX:age=3562.255340])129:18918.895955[&&NHX:age=8051.702368],((75:1545.314509[&&NHX:age=0.000000],24:1545.314509[&&NHX:age=0.000000])184:16499.310954[&&NHX:age=1545.314509],(74:12061.808515[&&NHX:age=0.000000],(34:2354.701987[&&NHX:age=0.000000],((64:1545.314509[&&NHX:age=0.000000],((81:0.000000[&&NHX:age=0.000000],97:0.000000[&&NHX:age=0.000000])194:0.000000[&&NHX:age=0.000000],79:0.000000[&&NHX:age=0.000000])162:1545.314509[&&NHX:age=0.000000])105:809.387478[&&NHX:age=1545.314509],(88:395.449492[&&NHX:age=0.000000],((((39:0.000000[&&NHX:age=0.000000],8:0.000000[&&NHX:age=0.000000])155:0.000000[&&NHX:age=0.000000],4:0.000000[&&NHX:age=0.000000])121:0.000000[&&NHX:age=0.000000],((89:0.000000[&&NHX:age=0.000000],84:0.000000[&&NHX:age=0.000000])178:0.000000[&&NHX:age=0.000000],2:0.000000[&&NHX:age=0.000000])168:0.000000[&&NHX:age=0.000000])113:0.000000[&&NHX:age=0.000000],(((20:0.000000[&&NHX:age=0.000000],((41:0.000000[&&NHX:age=0.000000],(71:0.000000[&&NHX:age=0.000000],(58:0.000000[&&NHX:age=0.000000],((80:0.000000[&&NHX:age=0.000000],(28:0.000000[&&NHX:age=0.000000],45:0.000000[&&NHX:age=0.000000])111:0.000000[&&NHX:age=0.000000])150:0.000000[&&NHX:age=0.000000],23:0.000000[&&NHX:age=0.000000])179:0.000000[&&NHX:age=0.000000])116:0.000000[&&NHX:age=0.000000])142:0.000000[&&NHX:age=0.000000])163:0.000000[&&NHX:age=0.000000],6:0.000000[&&NHX:age=0.000000])181:0.000000[&&NHX:age=0.000000])157:0.000000[&&NHX:age=0.000000],((49:0.000000[&&NHX:age=0.000000],44:0.000000[&&NHX:age=0.000000])195:0.000000[&&NHX:age=0.000000],((68:0.000000[&&NHX:age=0.000000],(66:0.000000[&&NHX:age=0.000000],50:0.000000[&&NHX:age=0.000000])132:0.000000[&&NHX:age=0.000000])136:0.000000[&&NHX:age=0.000000],(((17:0.000000[&&NHX:age=0.000000],47:0.000000[&&NHX:age=0.000000])104:0.000000[&&NHX:age=0.000000],12:0.000000[&&NHX:age=0.000000])187:0.000000[&&NHX:age=0.000000],((31:0.000000[&&NHX:age=0.000000],(25:0.000000[&&NHX:age=0.000000],16:0.000000[&&NHX:age=0.000000])197:0.000000[&&NHX:age=0.000000])106:0.000000[&&NHX:age=0.000000],(38:0.000000[&&NHX:age=0.000000],((82:0.000000[&&NHX:age=0.000000],((72:0.000000[&&NHX:age=0.000000],83:0.000000[&&NHX:age=0.000000])182:0.000000[&&NHX:age=0.000000],(59:0.000000[&&NHX:age=0.000000],93:0.000000[&&NHX:age=0.000000])118:0.000000[&&NHX:age=0.000000])166:0.000000[&&NHX:age=0.000000])164:0.000000[&&NHX:age=0.000000],(73:0.000000[&&NHX:age=0.000000],((69:0.000000[&&NHX:age=0.000000],48:0.000000[&&NHX:age=0.000000])138:0.000000[&&NHX:age=0.000000],5:0.000000[&&NHX:age=0.000000])191:0.000000[&&NHX:age=0.000000])146:0.000000[&&NHX:age=0.000000])153:0.000000[&&NHX:age=0.000000])151:0.000000[&&NHX:age=0.000000])125:0.000000[&&NHX:age=0.000000])185:0.000000[&&NHX:age=0.000000])117:0.000000[&&NHX:age=0.000000])175:0.000000[&&NHX:age=0.000000])177:0.000000[&&NHX:age=0.000000],(62:0.000000[&&NHX:age=0.000000],(76:0.000000[&&NHX:age=0.000000],(56:0.000000[&&NHX:age=0.000000],((46:0.000000[&&NHX:age=0.000000],(14:0.000000[&&NHX:age=0.000000],(9:0.000000[&&NHX:age=0.000000],(77:0.000000[&&NHX:age=0.000000],(15:0.000000[&&NHX:age=0.000000],(53:0.000000[&&NHX:age=0.000000],(((87:0.000000[&&NHX:age=0.000000],(95:0.000000[&&NHX:age=0.000000],(35:0.000000[&&NHX:age=0.000000],(10:0.000000[&&NHX:age=0.000000],91:0.000000[&&NHX:age=0.000000])186:0.000000[&&NHX:age=0.000000])100:0.000000[&&NHX:age=0.000000])190:0.000000[&&NHX:age=0.000000])174:0.000000[&&NHX:age=0.000000],(57:0.000000[&&NHX:age=0.000000],29:0.000000[&&NHX:age=0.000000])114:0.000000[&&NHX:age=0.000000])144:0.000000[&&NHX:age=0.000000],(26:0.000000[&&NHX:age=0.000000],(18:0.000000[&&NHX:age=0.000000],((65:0.000000[&&NHX:age=0.000000],(42:0.000000[&&NHX:age=0.000000],(98:0.000000[&&NHX:age=0.000000],61:0.000000[&&NHX:age=0.000000])122:0.000000[&&NHX:age=0.000000])167:0.000000[&&NHX:age=0.000000])130:0.000000[&&NHX:age=0.000000],3:0.000000[&&NHX:age=0.000000])135:0.000000[&&NHX:age=0.000000])141:0.000000[&&NHX:age=0.000000])103:0.000000[&&NHX:age=0.000000])115:0.000000[&&NHX:age=0.000000])123:0.000000[&&NHX:age=0.000000])139:0.000000[&&NHX:age=0.000000])154:0.000000[&&NHX:age=0.000000])137:0.000000[&&NHX:age=0.000000])109:0.000000[&&NHX:age=0.000000])183:0.000000[&&NHX:age=0.000000],0:0.000000[&&NHX:age=0.000000])131:0.000000[&&NHX:age=0.000000])112:0.000000[&&NHX:age=0.000000])152:0.000000[&&NHX:age=0.000000])124:0.000000[&&NHX:age=0.000000])145:0.000000[&&NHX:age=0.000000])161:395.449492[&&NHX:age=0.000000])193:1959.252495[&&NHX:age=395.449492])128:0.000000[&&NHX:age=2354.701987])196:9707.106528[&&NHX:age=2354.701987])176:5982.816948[&&NHX:age=12061.808515])189:8925.972860[&&NHX:age=18044.625462])148[&&NHX:age=26970.598323];";

//		input ="(((27:1002.805899[&&NHX:age=0.000000],(86:0.000000[&&NHX:age=0.000000],13:0.000000[&&NHX:age=0.000000])172:1002.805899[&&NHX:age=0.000000])107:0.000000[&&NHX:age=1002.805899],1:1002.805899[&&NHX:age=0.000000])101:25967.792423[&&NHX:age=1002.805899],(((75:1545.314509[&&NHX:age=0.000000],(96:1002.805899[&&NHX:age=0.000000],24:1002.805899[&&NHX:age=0.000000])192:542.508609[&&NHX:age=1002.805899])184:16499.310954[&&NHX:age=1545.314509],(22:12061.808515[&&NHX:age=0.000000],(92:12061.808515[&&NHX:age=0.000000],(74:12061.808515[&&NHX:age=0.000000],(52:8051.702368[&&NHX:age=0.000000],((88:395.449492[&&NHX:age=0.000000],4:395.449492[&&NHX:age=0.000000])193:1959.252495[&&NHX:age=395.449492],(64:122.586947[&&NHX:age=0.000000],((89:0.000000[&&NHX:age=0.000000],84:0.000000[&&NHX:age=0.000000])178:0.000000[&&NHX:age=0.000000],2:0.000000[&&NHX:age=0.000000])168:122.586947[&&NHX:age=0.000000])113:2232.115040[&&NHX:age=122.586947])128:5697.000381[&&NHX:age=2354.701987])119:4010.106147[&&NHX:age=8051.702368])176:0.000000[&&NHX:age=12061.808515])129:0.000000[&&NHX:age=12061.808515])173:5982.816948[&&NHX:age=12061.808515])189:8925.972860[&&NHX:age=18044.625462],((33:3562.255340[&&NHX:age=0.000000],(((81:0.000000[&&NHX:age=0.000000],97:0.000000[&&NHX:age=0.000000])194:0.000000[&&NHX:age=0.000000],79:0.000000[&&NHX:age=0.000000])162:395.449492[&&NHX:age=0.000000],(19:49.193481[&&NHX:age=0.000000],((17:0.000000[&&NHX:age=0.000000],47:0.000000[&&NHX:age=0.000000])104:0.000000[&&NHX:age=0.000000],12:0.000000[&&NHX:age=0.000000])187:49.193481[&&NHX:age=0.000000])185:346.256011[&&NHX:age=49.193481])149:3166.805848[&&NHX:age=395.449492])105:1801.590881[&&NHX:age=3562.255340],((54:3562.255340[&&NHX:age=0.000000],(34:49.193481[&&NHX:age=0.000000],(39:0.000000[&&NHX:age=0.000000],8:0.000000[&&NHX:age=0.000000])155:49.193481[&&NHX:age=0.000000])121:3513.061859[&&NHX:age=49.193481])108:1801.590881[&&NHX:age=3562.255340],((((94:0.000000[&&NHX:age=0.000000],((90:0.000000[&&NHX:age=0.000000],67:0.000000[&&NHX:age=0.000000])180:0.000000[&&NHX:age=0.000000],(51:0.000000[&&NHX:age=0.000000],40:0.000000[&&NHX:age=0.000000])102:0.000000[&&NHX:age=0.000000])134:0.000000[&&NHX:age=0.000000])188:0.000000[&&NHX:age=0.000000],(43:0.000000[&&NHX:age=0.000000],((((99:0.000000[&&NHX:age=0.000000],70:0.000000[&&NHX:age=0.000000])198:0.000000[&&NHX:age=0.000000],37:0.000000[&&NHX:age=0.000000])140:0.000000[&&NHX:age=0.000000],(85:0.000000[&&NHX:age=0.000000],55:0.000000[&&NHX:age=0.000000])110:0.000000[&&NHX:age=0.000000])170:0.000000[&&NHX:age=0.000000],36:0.000000[&&NHX:age=0.000000])147:0.000000[&&NHX:age=0.000000])171:0.000000[&&NHX:age=0.000000])159:0.000000[&&NHX:age=0.000000],((60:0.000000[&&NHX:age=0.000000],(63:0.000000[&&NHX:age=0.000000],(78:0.000000[&&NHX:age=0.000000],(30:0.000000[&&NHX:age=0.000000],21:0.000000[&&NHX:age=0.000000])165:0.000000[&&NHX:age=0.000000])156:0.000000[&&NHX:age=0.000000])126:0.000000[&&NHX:age=0.000000])120:0.000000[&&NHX:age=0.000000],7:0.000000[&&NHX:age=0.000000])169:0.000000[&&NHX:age=0.000000])143:1545.314509[&&NHX:age=0.000000],(11:639.178343[&&NHX:age=0.000000],(((20:0.000000[&&NHX:age=0.000000],((41:0.000000[&&NHX:age=0.000000],(71:0.000000[&&NHX:age=0.000000],(58:0.000000[&&NHX:age=0.000000],((80:0.000000[&&NHX:age=0.000000],(28:0.000000[&&NHX:age=0.000000],45:0.000000[&&NHX:age=0.000000])111:0.000000[&&NHX:age=0.000000])150:0.000000[&&NHX:age=0.000000],23:0.000000[&&NHX:age=0.000000])179:0.000000[&&NHX:age=0.000000])116:0.000000[&&NHX:age=0.000000])142:0.000000[&&NHX:age=0.000000])163:0.000000[&&NHX:age=0.000000],6:0.000000[&&NHX:age=0.000000])181:0.000000[&&NHX:age=0.000000])157:0.000000[&&NHX:age=0.000000],((49:0.000000[&&NHX:age=0.000000],44:0.000000[&&NHX:age=0.000000])195:0.000000[&&NHX:age=0.000000],((68:0.000000[&&NHX:age=0.000000],(66:0.000000[&&NHX:age=0.000000],50:0.000000[&&NHX:age=0.000000])132:0.000000[&&NHX:age=0.000000])136:0.000000[&&NHX:age=0.000000],((31:0.000000[&&NHX:age=0.000000],(25:0.000000[&&NHX:age=0.000000],16:0.000000[&&NHX:age=0.000000])197:0.000000[&&NHX:age=0.000000])106:0.000000[&&NHX:age=0.000000],(38:0.000000[&&NHX:age=0.000000],((82:0.000000[&&NHX:age=0.000000],((72:0.000000[&&NHX:age=0.000000],83:0.000000[&&NHX:age=0.000000])182:0.000000[&&NHX:age=0.000000],(59:0.000000[&&NHX:age=0.000000],93:0.000000[&&NHX:age=0.000000])118:0.000000[&&NHX:age=0.000000])166:0.000000[&&NHX:age=0.000000])164:0.000000[&&NHX:age=0.000000],(73:0.000000[&&NHX:age=0.000000],((69:0.000000[&&NHX:age=0.000000],48:0.000000[&&NHX:age=0.000000])138:0.000000[&&NHX:age=0.000000],5:0.000000[&&NHX:age=0.000000])191:0.000000[&&NHX:age=0.000000])146:0.000000[&&NHX:age=0.000000])153:0.000000[&&NHX:age=0.000000])151:0.000000[&&NHX:age=0.000000])125:0.000000[&&NHX:age=0.000000])117:0.000000[&&NHX:age=0.000000])175:0.000000[&&NHX:age=0.000000])177:0.000000[&&NHX:age=0.000000],(62:0.000000[&&NHX:age=0.000000],(76:0.000000[&&NHX:age=0.000000],(56:0.000000[&&NHX:age=0.000000],((46:0.000000[&&NHX:age=0.000000],(14:0.000000[&&NHX:age=0.000000],(9:0.000000[&&NHX:age=0.000000],(77:0.000000[&&NHX:age=0.000000],(15:0.000000[&&NHX:age=0.000000],(53:0.000000[&&NHX:age=0.000000],(((87:0.000000[&&NHX:age=0.000000],(95:0.000000[&&NHX:age=0.000000],(35:0.000000[&&NHX:age=0.000000],(10:0.000000[&&NHX:age=0.000000],91:0.000000[&&NHX:age=0.000000])186:0.000000[&&NHX:age=0.000000])100:0.000000[&&NHX:age=0.000000])190:0.000000[&&NHX:age=0.000000])174:0.000000[&&NHX:age=0.000000],(57:0.000000[&&NHX:age=0.000000],29:0.000000[&&NHX:age=0.000000])114:0.000000[&&NHX:age=0.000000])144:0.000000[&&NHX:age=0.000000],(26:0.000000[&&NHX:age=0.000000],((32:0.000000[&&NHX:age=0.000000],18:0.000000[&&NHX:age=0.000000])127:0.000000[&&NHX:age=0.000000],((65:0.000000[&&NHX:age=0.000000],(42:0.000000[&&NHX:age=0.000000],(98:0.000000[&&NHX:age=0.000000],61:0.000000[&&NHX:age=0.000000])122:0.000000[&&NHX:age=0.000000])167:0.000000[&&NHX:age=0.000000])130:0.000000[&&NHX:age=0.000000],3:0.000000[&&NHX:age=0.000000])135:0.000000[&&NHX:age=0.000000])141:0.000000[&&NHX:age=0.000000])103:0.000000[&&NHX:age=0.000000])115:0.000000[&&NHX:age=0.000000])123:0.000000[&&NHX:age=0.000000])139:0.000000[&&NHX:age=0.000000])154:0.000000[&&NHX:age=0.000000])137:0.000000[&&NHX:age=0.000000])109:0.000000[&&NHX:age=0.000000])183:0.000000[&&NHX:age=0.000000],0:0.000000[&&NHX:age=0.000000])131:0.000000[&&NHX:age=0.000000])112:0.000000[&&NHX:age=0.000000])152:0.000000[&&NHX:age=0.000000])124:0.000000[&&NHX:age=0.000000])145:639.178343[&&NHX:age=0.000000])161:906.136166[&&NHX:age=639.178343])133:3818.531712[&&NHX:age=1545.314509])196:0.000000[&&NHX:age=5363.846221])160:21606.752102[&&NHX:age=5363.846221])148:0.000000[&&NHX:age=26970.598323])158[&&NHX:age=26970.598323];";

//		input ="(((4:5363.846221[&&NHX:age=0.000000],7:5363.846221[&&NHX:age=0.000000])13:2687.856147[&&NHX:age=5363.846221],(6:3562.255340[&&NHX:age=0.000000],1:3562.255340[&&NHX:age=0.000000])14:4489.447028[&&NHX:age=3562.255340])12:0.000000[&&NHX:age=8051.702368],(3:5363.846221[&&NHX:age=0.000000],(0:2354.701987[&&NHX:age=0.000000],(2:2354.701987[&&NHX:age=0.000000],5:2354.701987[&&NHX:age=0.000000])9:0.000000[&&NHX:age=2354.701987])10:3009.144234[&&NHX:age=2354.701987])11:2687.856147[&&NHX:age=5363.846221])8[&&NHX:age=8051.702368];";
		String[] inputArray = input.split("\\(|,|\\)");

		ArrayList<String> temp = new ArrayList<String>();

		for(int i = 0; i < inputArray.length; i++){
			if(inputArray[i].length() > 0){
				temp.add(inputArray[i]);
			}
		}

		//		System.out.println("" + temp.size());

		Double[][] age = new Double[temp.size() - 1][4];

		for(int i = 0; i < temp.size() - 1; i++){
			String[] currentThing = temp.get(i).split(":|\\[");
			String[] ageArray = currentThing[3].split("=|\\]");
			age[i][1] = Double.parseDouble(currentThing[0]); // id
			age[i][0] = Double.parseDouble(ageArray[1]); // age
			age[i][2] = Double.parseDouble(currentThing[1]); // branch length
		}

		Arrays.sort(age, new Comparator<Double[]>() {
			@Override
			public int compare(final Double[] entry1, final Double[] entry2){

				int value = entry1[0].compareTo(entry2[0]);

				if (value == 0){
					value = entry1[2].compareTo(entry2[2]);
				}

				if(value == 0){
					value = entry1[1].compareTo(entry2[1]);
				}

				return value;
			}

		});

		String[] mcraArray = temp.get(temp.size()-1).split("\\[");
		int rootBranch = Integer.parseInt(mcraArray[0]);

		int index = input.indexOf(")" + rootBranch + "[");

		Node root = fillNode(new Node(), rootBranch, index, input, age);

		Integer numTaxa = (temp.size() + 1) / 2;

		int[][] tree = new int[(temp.size() - numTaxa.intValue())][3];

		for(int i = tree.length - 1; i >= 0; i--){
//			System.out.println(i);
			ArrayList<Node> nodeList = this.nextCoalese(root, numTaxa);
			Collections.shuffle(nodeList);
			
			double maxAge = -2.0;
			int ageIndex = -1;
			
			for(int j = 0; j < nodeList.size(); j++){
				if (nodeList.get(j).age > maxAge){
					maxAge = nodeList.get(j).age;
					ageIndex = j;
				}
			}
			
			Node currentNode = nodeList.get(ageIndex);
			tree[i][2] = currentNode.key;
			tree[i][0] = currentNode.lChild;
			tree[i][1] = currentNode.rChild;
			currentNode.coalesed = true;
		}

		for(int i = 0; i < tree.length; i++){
			System.out.println("" + tree[i][0] + " " + tree[i][1] + " " + tree[i][2]);
		}
		
		System.out.println("---");
		
		for(int i = 0; i < tree.length; i++){
			int currentBranch = tree[i][2];
			tree[i][2] = i + numTaxa;
			for(int j = i + 1; j < tree.length; j++){
				if(tree[j][0] == currentBranch){
					tree[j][0] = i + numTaxa;
					break;
				}
				if(tree[j][1] == currentBranch){
					tree[j][1] = i + numTaxa;
					break;
				}
			}
		}
		
		return tree;
	}

	public Node fillNode(Node currentNode, int currentBranch, int currentIndex, String input,Double age[][]){
		currentNode.key = currentBranch;

		int[] kids = findChildren(input, currentIndex);

		currentNode.lChild = kids[0];
		currentNode.rChild = kids[1];

		for(int i = 0; i < age.length; i++){
			if (age[i][1].intValue() == currentBranch){
//				System.out.println("Found age for " + currentBranch + " " + age[i][0]);
				currentNode.age = age[i][0];
			}
		}

		int child1Index = findIndex(input, kids[0]);
		int child2Index = findIndex(input, kids[1]);

		if(kids[0] > -1){
			currentNode.leftChild = fillNode(new Node(), kids[0], child1Index, input, age);
		}

		if(kids[1] > -1){
			currentNode.rightChild = fillNode(new Node(), kids[1], child2Index, input, age);
		}

		return currentNode;
	}

	public int[] findChildren(String input, int index){

		int parens = 0;
		int found = 0;

		int[] children = new int[2];
		children[0] = -1;
		children[1] = -1;

		for (int j = index; j >= 0; j--){
			if (input.charAt(j) == ')'){
				parens++;
			}
			if (input.charAt(j) == '('){
				parens--;
			}
			if(parens == 1 && input.charAt(j) == ':' && input.charAt(j + 1) != 'a'){

				//				System.out.println(input.substring(j - 2, j));
				if(j > 3 && input.substring(j - 4, j).matches("\\d{4}")){
					//							System.out.println(input.substring(j - 2, j));
					if(found == 0){
						children[0] = Integer.parseInt(input.substring(j - 4, j));
					} else if(found == 1){
						children[1] = Integer.parseInt(input.substring(j - 4, j));
					}
					found++;

				} else if(j > 2 && input.substring(j - 3, j).matches("\\d{3}")){
					//							System.out.println(input.substring(j - 2, j));
					if(found == 0){
						children[0] = Integer.parseInt(input.substring(j - 3, j));
					} else if(found == 1){
						children[1] = Integer.parseInt(input.substring(j - 3, j));
					}
					found++;

				} else if(j > 1 && input.substring(j - 2, j).matches("\\d{2}")){
					//							System.out.println(input.substring(j - 2, j));
					if(found == 0){
						children[0] = Integer.parseInt(input.substring(j - 2, j));
					} else if(found == 1){
						children[1] = Integer.parseInt(input.substring(j - 2, j));
					}
					found++;

				} else if(j > 0 && input.substring(j - 1, j).matches("\\d{1}")){
					if(found == 0){
						children[0] = Integer.parseInt(input.substring(j - 1, j));
					} else if(found == 1){
						children[1] = Integer.parseInt(input.substring(j - 1, j));
					}
					found++;
				}
				

			}

			if(found > 1){
				break;
			}
		}

		return children;

	}

	public int findIndex(String input, int currentBranch){
		int result = input.indexOf(")" + currentBranch + ":");
		

		return result;
	}

	public ArrayList<Node> nextCoalese(Node currentNode, int numTaxa){
		if(currentNode.key < numTaxa){
			ArrayList<Node> returnList = new ArrayList<Node>();
			return returnList;
		}
		
		if(currentNode.coalesed == false){
			ArrayList<Node> returnList = new ArrayList<Node>();
			returnList.add(currentNode);
			return returnList;
		}

		ArrayList<Node> leftChild = new ArrayList<Node>();
		if(currentNode.leftChild != null){
			ArrayList<Node> temp = nextCoalese(currentNode.leftChild, numTaxa);
			leftChild.addAll(temp);
		}
		ArrayList<Node> rightChild = new ArrayList<Node>();
		if(currentNode.rightChild != null){
			ArrayList<Node> temp = nextCoalese(currentNode.rightChild, numTaxa);
			rightChild.addAll(temp);
		}
		
		leftChild.addAll(rightChild);

		return leftChild;

	}
	
	public ArrayList<Integer> getMembers(int[][] tree, int branch) {
		ArrayList<Integer> members = new ArrayList<Integer>();

		if (tree[branch][0] < this.numberOfTaxa) {
			Integer member1 = tree[branch][0];
			members.add(member1);
		} else {
			members.addAll(this.getMembers(tree, tree[branch][0] - this.numberOfTaxa));
		}

		if (tree[branch][1] < this.numberOfTaxa) {
			Integer member2 = tree[branch][1];
			members.add(member2);
		} else {
			members.addAll(this.getMembers(tree, tree[branch][1] - this.numberOfTaxa));
		}

		return members;
	}

}