#include <stdio.h>

typedef int dtype;

int main(int argc, char *argv[])
{
	int n;
	int j;	    		//ѭ�����Ʊ���
	FILE *fp;
	float a;
	printf("Please type in the number of values:");
	scanf("%d ",&n);
	if((fp = fopen(argv[1],"wb")) == NULL)
	{
		printf("Cannot open file!");
		//���ﻹ��֪�����ʲô�˳���� exit(1);
	}				//��ָ���ļ�
				
	for(j = 0; j < n; j++)
	{
		scanf("%f", &a);
		fwrite(&a,sizeof(int),1,fp);
	}
	return 0;	
}
