1.使用OpenBabel将pdb进行转换

obabel $PDB -O $MOL2
obabel $MOL2 -O $PDB

2.利用list_gen.py生成list

3.使用make_ploy_lt生成lt文件

4.写system.lt文件

SiO paper：A Force Field and a Surface Model Database for Silica to Simulate Interfacial Properties in Atomic Resolution


moltemplate.sh -pdb box_ob.pdb -atomstyle full system.lt