#include <stdio.h>

typedef int dtype;

int main(int argc, char *argv[])
{
	int n;
	int j;	    		//循环控制变量
	FILE *fp;
	float a;
	printf("Please type in the number of values:");
	scanf("%d ",&n);
	if((fp = fopen(argv[1],"wb")) == NULL)
	{
		printf("Cannot open file!");
		//这里还不知道添加什么退出语句 exit(1);
	}				//打开指定文件
				
	for(j = 0; j < n; j++)
	{
		scanf("%f", &a);
		fwrite(&a,sizeof(int),1,fp);
	}
	return 0;	
}
