#include<iostream>
#include<graphics.h>
#include<cstdlib>
#include<cstring>
#include<algorithm>
#include<conio.h>
#include<ctime>
#include<cmath>
#include<cstdio>
#include<stack>
#include<string>

using namespace std;

/*设计思想：
随机选取点添边，同时构建并查集检验图的连通性，若生成的图非连通则添边，若仍不满足条件，则重新构建直到满足条件。
通过Floyd算法获得距离数组和路径数组，用于后续匹配和查找。获得奇点下标以及奇点个数，用状压dp进行最小权匹配，获得奇点对。
在原图中添加边消除奇点获得欧拉图，用DFS找到一条欧拉回路。在欧拉回路中，完善不直接相连奇点的路径，获得实际路径。
获取点的坐标，固定邮局A点，根据结点个数确定权值，判断两点间距离，如果太小则重新选取，最终得到分布均匀的图。
画出所有结点，令点按路径移动，最终连线，动态演示完毕。
*/

#define INF 65535			//定义无穷大
#define MaxVerNum  100		//最大结点数目
#define eps 1e-3			//精度，在精度之内则符合要求

string num_to_uppercase(int num) {
    // 定义数字和字母的对应关系
    char num_to_letter[] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};
    // 检查数字是否在 1-26 范围内，否则返回空字符串
    if (num < 1 || num > 26) {
        return "";
    }
    // 将数字转换为大写字母
    return string(1, num_to_letter[num - 1]);
}

typedef char elementType;	//定义图中顶点的数据类型
typedef int eInfoType;		//定义 eInfo 的数据类型，即权值的数据类型

typedef struct eNode		//定义边链表的结点结构
{
	int adjVer;				//邻接顶点信息，此处为顶点编号，从 1 开始
	eInfoType eInfo;		//边链表中表示边的相关信息，比如表的权值
	struct eNode* next;		//指向边链表中的下一个结点。
	bool visited;
}EdgeNode;					//边链表结点类型

typedef struct vNode		//定义顶点表的结点结构
{
	elementType data;		//顶点表数据类型
	EdgeNode* firstEdge;	//指向边结点
}VerNode;					//顶点表结点类型

typedef struct oNode {		//奇点对存储结构
	int node_1;
	int node_2;
}OddNode;

typedef struct Euler {		//路径结点定义
	int point;				//结点下标
	struct Euler* next;		//指向下一个路径结点
}Node;

typedef  struct {			//定义结点基本信息的结构
	int x, y;				//横、纵坐标
	elementType name;		//结点名称A ~ Z
	bool exist;				//是否已获得坐标
}vN;

stack<int> Q;													//全局栈，用于获得结点之间的路径
int bg, to;														//边的起点和终点
int VerNum;														//顶点数
int odd_num = 0, lf = 0;										//奇点个数
int* dp;														//全局指针，根据奇点个数动态分派内存
bool isIP = true;												//用于判断VerNum是否已有值
int degree[MaxVerNum];											//表示结点度数
int tree[MaxVerNum];											//并查集，用于判断图的连通性
int path[MaxVerNum][MaxVerNum];									//路径数组
VerNode VerList[MaxVerNum];										//顶点表
eInfoType dist[MaxVerNum][MaxVerNum];							//距离矩阵
eInfoType AdjMatrix[MaxVerNum][MaxVerNum];						//邻接矩阵
int* odd_point;													//记录奇点
int odd_n[MaxVerNum];											//用于记录奇点个数及奇点是第几个结点
OddNode odd_match[MaxVerNum];									//奇点间配对
void initialize();												//初始化函数
int findroot(int x);											//建立并查集
bool add(int m, int n, int w, bool instant);					//添加边
void free();													//释放顶点表的边链表及动态创建的奇点表
void connected();												//若生成的图非连通，则合并
bool Legal();													//判断图是否符合条件
void odd();														//不符合条件则重新生成图
void Floyd();													//获得任意两点距离矩阵以及路径
int get_oddedge(int x);											//获得需要添边的奇点对
void extend_edge();												//在原图中添加边以获得欧拉图
bool IsExtend(int x, int y);									//判断是否是添加的边
void search_euler(int u);										//寻找一条欧拉回路并记录
void get_path(Node* p);											//获得实际路径
void getarc(int& arc, int x, int y);							//根据点的位置获得合适的角度
void getWeight(int& weight);									//根据点的数目确定点与点的距离
void getXY(vN* xy);												//获取点的坐标
void generate(vN* xy, int i, EdgeNode* p);						//产生点的坐标
void adjustXY(vN* xy);											//调整点的坐标，使各点之间保持一定距离
void paint(vN* xy);												//绘制出所有结点
void drawline(vN* xy);											//绘制点之间的边
void show(vN* xy, Node* Euler);									//动态演示求解过程（路径）
void getd(double& d, double k);									//因斜率不同选取合适的步长
void flag(int x, int y);										//到达一个结点则高亮一次，用于确定邮递员的位置
void move(int x1, int y1, int x2, int y2, vN* xy);				//演示邮递员在图上的移动过程

void initialize()																	//初始化函数
{
	int count = 0, i, m, n, weight, EdgeNum;
	bool succeed;
	if (isIP)																		//是否需要再次生成结点个数
	{
		VerNum = rand() % 16 + 9;													//结点个数在9 ~ 24
		isIP = false;
	}
	for (i = 1; i < VerNum + 1; i++)												//初始化顶点表
	{
		VerList[i].data = char(i - 1) + 'A';
		VerList[i].firstEdge = NULL;
	}
	odd_point = (int*)malloc((VerNum + 1) * sizeof(int));						//动态数组
	memset(degree, 0, sizeof(degree));									//初始化dgree数组
	memset(tree, -1, sizeof(tree));									//初始化tree数组
	EdgeNum = rand() % (VerNum * 2 - VerNum + 1) + VerNum - 1;						//控制边的个数
	while (count < EdgeNum)															//确定两结点以及权值
	{
		m = rand() % VerNum + 1;
		n = rand() % VerNum + 1;
		weight = rand() % 15 + 4;
		while (m == n)
		{
			n = rand() % VerNum + 1;
		}
		succeed = add(m, n, weight, false);
		if (succeed)
			count++;
	}
}

int findroot(int x)																	//并查集的方式判断图连通性
{
	if (tree[x] == -1)
		return x;
	else
	{
		int temp = findroot(tree[x]);
		tree[x] = temp;
		return temp;
	}
}

bool add(int m, int n, int weight, bool instant)									//添边，根据instant的值来确定何种形式的添边
{
	EdgeNode* p = NULL, * edge1 = NULL, * edge2 = NULL;
	int t, tempa, tempb;
	if (!instant)																	//不允许重复边的添边方式
	{																				//不符合条件则返回false			
		if (degree[m] > 3 && degree[n] > 3 && m < 1 || m > VerNum || n < 1 || n > VerNum)
			return false;
		p = VerList[m].firstEdge;													//防重
		if (p)
		{
			while (p->next && p->adjVer != n)
			{
				p = p->next;
			}
		}
		if (p)
		{
			if (p->adjVer == n)
			{
				while (1)
				{
					t = n;
					p = VerList[m].firstEdge;
					while (n == m || n == t)
						n = rand() % VerNum + 1;
					while (p->next && p->adjVer != n)
					{
						p = p->next;
					}
					if (p->adjVer != n)
						break;
				}
			}
		}
		edge1 = new EdgeNode;														//添加双向边
		edge1->adjVer = n;
		edge1->eInfo = weight;
		edge1->next = NULL;
		edge1->visited = false;
		if (VerList[m].firstEdge == NULL)
			VerList[m].firstEdge = edge1;
		else
			p->next = edge1;
		p = VerList[n].firstEdge;
		if (p)
		{
			while (p->next)
			{
				p = p->next;
			}
		}
		edge2 = new EdgeNode;
		edge2->adjVer = m;
		edge2->eInfo = weight;
		edge2->next = NULL;
		edge2->visited = false;
		if (VerList[n].firstEdge == NULL)
			VerList[n].firstEdge = edge2;
		else
			p->next = edge2;
		tempa = findroot(m);														//建立并查集用以判断连通性
		tempb = findroot(n);
		if (tempa != tempb)															//若祖先不同则合并子图
			tree[tempa] = tempb;
	}
	else                                                                            //重复边的添边
	{
		p = VerList[m].firstEdge;
		while (p->next)
		{
			p = p->next;
		}
		edge1 = new EdgeNode;
		edge1->adjVer = n;
		edge1->eInfo = weight;
		edge1->next = NULL;
		edge1->visited = false;
		p->next = edge1;

		p = VerList[n].firstEdge;
		while (p->next)
		{
			p = p->next;
		}
		edge2 = new EdgeNode;
		edge2->adjVer = m;
		edge2->eInfo = weight;
		edge2->next = NULL;
		edge2->visited = false;
		p->next = edge2;
	}
	degree[m]++;																	//结点度数相应增加
	degree[n]++;
	return true;
}

void connected()
{
	int i, ans = 0, a[80], weight;
	memset(a, -1, sizeof(a));											//初始化用以存放独立子图祖先结点的数组
	for (i = 1; i < VerNum + 1; i++)
	{																				//判断连通性
		if (tree[i] == -1)
		{
			a[ans] = i;
			ans++;
		}
	}
	i = 0;
	if (ans > 1)																	//合并子图，最后仅有一个连通图
	{
		while (i < ans - 1)
		{
			weight = rand() % 15 + 4;
			add(a[i], a[i + 1], weight, false);
			i++;
		}
	}
}

bool Legal()																		//判断子图是否合法
{
	int i, odd = 0;
	for (i = 1; i < VerNum + 1; i++)												//统计奇点个数
	{
		if (degree[i] % 2)
			odd++;
	}

	if (odd == 0 || odd > 2)
		return true;
	if (odd == 2)
	{
		if (degree[1] % 2 == 1)														//当奇点个数为2，若1不为奇点，则不符合要求
			return true;
		else
			return false;
	}
}

void odd()																			//循环构建图，直到生成的图符合要求
{
	while (!Legal())
	{
		free();
		initialize();
		connected();
	}
}

void Floyd()																		//Floyd算法确定任意两点距离及路径
{
	int i, j, k;
	EdgeNode* p = NULL;
	for (i = 1; i < VerNum + 1; i++)
		for (j = 1; j < VerNum + 1; j++)
			AdjMatrix[i][j] = INF;
	for (i = 1; i < VerNum + 1; i++)												//将邻接表转换为邻接矩阵
	{
		p = VerList[i].firstEdge;
		while (p)
		{
			AdjMatrix[i][p->adjVer] = p->eInfo;
			p = p->next;
		}
	}
	for (i = 1; i < VerNum + 1; i++)												//初始化路径数组
	{
		for (j = 1; j < VerNum + 1; j++)
		{
			dist[i][j] = AdjMatrix[i][j];
			if (i != j && AdjMatrix[i][j] < INF)
				path[i][j] = i;
			else
				path[i][j] = -1;
		}
	}
	for (k = 1; k < VerNum + 1; k++)
		for (i = 1; i < VerNum + 1; i++)
			for (j = 1; j < VerNum + 1; j++)
			{
				if (i != j && dist[i][k] + dist[k][j] < dist[i][j])
				{
					dist[i][j] = dist[i][k] + dist[k][j];
					path[i][j] = path[k][j];
				}
			}
}

int get_oddedge(int x)																//获得能最小权匹配的奇点对
{
	int i, j;
	int ans = INF, temp, ANS;
	EdgeNode* point = NULL;
	if (dp[x] || !x)
		return dp[x];
	for (i = 0; i < odd_num - 1; i++)												//遍历压缩状态中的所有'1'对，即所有的奇点对
	{
		if ((1 << i) & x)
		{
			for (j = i + 1; j < odd_num; j++)
			{
				if ((1 << j) & x)
				{
					temp = x ^ (1 << i) ^ (1 << j);
					ANS = get_oddedge(temp) + dist[odd_n[odd_num - i]][odd_n[odd_num - j]];
					if (ans > ANS)
					{
						ans = ANS;													//记录权值最小的奇点对下标
						bg = odd_n[odd_num - i];
						to = odd_n[odd_num - j];
					}
				}
			}
		}
	}
	return dp[x] = ans;
}

void extend_edge()																	//添加奇点对的边
{
	int i, j = 0, a = 0, count = 1, n = 0;
	bool instant = false;
	EdgeNode* point = NULL;
	cout << "奇点有：" << endl;
	for (i = 1; i < VerNum + 1; i++)												//初始化奇点数组
	{
		if (degree[i] % 2)
		{
			odd_num++;
			odd_n[odd_num] = i;
			odd_point[odd_num] = i;
			cout << num_to_uppercase(i) << "\t";
		}
	}
	cout << endl;
	lf = odd_num;
	cout << "奇点个数为：" << lf << endl;
	for (i = odd_num; i >= 1; i--)													//根据奇点个数生成对应二进制数
	{
		a <<= 1;
		a++;
	}
	i = 0;
	while (a)
	{
		dp = (int*)malloc((1 << odd_num) * sizeof(int));							//创建并初始化dp数组
		memset(dp, 0, (1 << odd_num) * sizeof(int));
		instant = false;
		get_oddedge(a);																//获取奇点对并记录
		odd_match[i].node_1 = bg;
		odd_match[i].node_2 = to;
		i++;
		for (j = 1; j <= odd_num; j++)												//更新奇点数组（添边后的奇点删除）
		{
			if (odd_n[j] != bg && odd_n[j] != to)
			{
				odd_n[count++] = odd_n[j];
			}
		}
		odd_num -= 2;
		count = 1;
		a = 0;
		for (int k = odd_num; k >= 1; k--)
		{
			a <<= 1;
			a++;
		}
		point = VerList[bg].firstEdge;												//添边
		while (point)																//判断两奇点是否直接相连
		{
			if (point->adjVer == to)
			{
				instant = true;
				break;
			}
			point = point->next;
		}
		add(bg, to, dist[bg][to], instant);
		//cout << bg << "\t" << to << "\t" << dist[bg][to] << endl;
		free(dp);																	//释放dp数组
	}
	for (int i = 1; i <= VerNum; i++)
	{
		if (degree[i] % 2)
			n++;
	}
	cout << "添边后的奇点个数为：" << n << endl;
}

bool IsExtend(int x, int y)															//是否为后来添加的边
{
	int i;
	for (i = 0; i < lf / 2; i++)
	{
		if (odd_match[i].node_1 == x || odd_match[i].node_2 == x)
		{
			if (odd_match[i].node_1 == y || odd_match[i].node_2 == y)
				return true;
		}
	}
	return false;
}

void search_euler(int u)															 //DFS寻找欧拉回路
{
	EdgeNode* p = VerList[u].firstEdge;;
	EdgeNode* q = NULL;
	while (p)
	{
		if (!p->visited)
		{
			p->visited = true;
			q = VerList[p->adjVer].firstEdge;
			while (q && q->adjVer != u)
			{
				q = q->next;
			}
			while (q && q->visited)													//防止重边无法正确遍历
			{
				q = q->next;
				while (q && q->adjVer != u)
				{
					q = q->next;
				}
			}
			if (!q)
				cout << "bug" << endl;
			if (q)
				q->visited = true;
			search_euler(p->adjVer);
		}
		p = p->next;
	}
	Q.push(u);																		//邻接点都被遍历则入栈
}

void get_path(Node* p)																//获得可行路径
{
	Node* q = p, * t = NULL, * r = NULL;
	int x, k;
	bool instant = true;
	while (!Q.empty())																//生成欧拉回路
	{
		t = new Node;
		x = Q.top();
		Q.pop();
		t->point = x;
		t->next = NULL;
		q->next = t;
		q = t;
	}
	if (p->next)
	{
		r = q = p->next;
		t = q->next;
	}
	else
	{
		cout << "此题无解" << endl;
		exit(1);
	}
	cout << "欧拉回路为：" << endl;
	while (r)
	{
		cout << char(r->point - 1 + 'A');
		if (r->next)
			cout << "->";
		r = r->next;
	}
	cout << endl;
	while (t)																	//得到最终路径
	{
		instant = true;
		if (IsExtend(q->point, t->point))										//补全非直接相连奇点间的路径
		{
			k = path[q->point][t->point];
			while (k != q->point)
			{
				Q.push(k);
				k = path[q->point][k];
				instant = false;
			}
			r = q;
			while (!Q.empty())
			{
				x = Q.top();
				Node* node = new Node;
				node->point = x;
				r->next = node;
				r = node;
				r->next = t;
				Q.pop();
			}
			r = NULL;
			if (!instant)
				t = q->next;
			else
			{
				q = q->next;
				t = q->next;
			}
		}
		else
		{
			q = t;
			t = t->next;
		}
	}
	if (p->next)
		q = p->next;
	cout << "最终路径为：" << endl;
	while (q)
	{
		cout << char(q->point - 1 + 'A');
		if (q->next)
			cout << "->";
		q = q->next;
	}
	cout << endl;
}

void free()																			//释放动态分配的内存
{
	int i;
	EdgeNode* p = NULL, * u = NULL;
	free(odd_point);
	for (i = 1; i < VerNum + 1; i++)
	{
		p = VerList[i].firstEdge;
		VerList[i].firstEdge = NULL;
		while (p)
		{
			u = p->next;
			delete p;
			p = u;
		}
	}
}

void getarc(int& arc, int x, int y)													//分区域获得角度，使分布均匀
{
	int n;
	if (x >= 75 && x <= 530)
	{
		if (y >= 45 && y < 288)
		{
			arc = rand() % 90 + 270;
		}
		else if (y >= 288 && y < 572)
		{
			n = rand() % 2 + 1;
			if (n == 1)
				arc = rand() % 90;
			else
				arc = rand() % 90 + 270;
		}
		else if (y >= 572 && y <= 815)
		{
			arc = rand() % 90;
		}
	}
	else if (x > 530 && x <= 985)
	{
		if (y >= 45 && y < 288)
		{
			arc = rand() % 180 + 180;
		}
		else if (y >= 288 && y < 572)
		{
			arc = rand() % 360;
		}
		else if (y >= 572 && y <= 815)
		{
			arc = rand() % 180;
		}
	}
	else if (x > 985 && x <= 1450)
	{
		if (y >= 45 && y <= 288)
		{
			arc = rand() % 90 + 180;
		}
		else if (y > 288 && y <= 572)
		{
			arc = rand() % 180 + 90;
		}
		else if (y > 572 && y <= 815)
		{
			arc = rand() % 90 + 90;
		}
	}
}

void getWeight(int& weight)															 //根据点的个数确定距离
{
	switch (VerNum)
	{
	case 24:
	case 23:
		weight = 140;
		break;
	case 22:
		weight = 145;
		break;
	case 21:
		weight = 150;
		break;
	case 20:
	case 19:
		weight = 155;
		break;
	case 18:
		weight = 160;
		break;
	case 17:
	case 16:;
	case 15:
		weight = 170;
		break;
	case 14:
	case 13:
	case 12:
	case 11:
	case 10:
	case 9:
		weight = 180;
		break;
	}
}

void generate(vN* xy, int i, EdgeNode* p)											//在以某点为圆心的圆上取得一点
{
	int weight, arc;
	int X = 0, Y = 0;
	if (!xy[p->adjVer].exist)
	{
		getWeight(weight);
		while (X < 75 || X > 1450 || Y < 45 || Y > 815)
		{
			getarc(arc, xy[i].x, xy[i].y);
			X = int(xy[i].x + weight * cos(arc * 3.14 / 180) + 0.5);
			Y = int(xy[i].y + weight * sin(arc * 3.14 / 180) + 0.5);
		}
		xy[p->adjVer].x = X;
		xy[p->adjVer].y = Y;
		xy[p->adjVer].exist = true;
	}
}

void getXY(vN* xy)																	//获得结点XY坐标
{
	int i, k;
	int X = 0, Y = 0;
	bool all;
	EdgeNode* p = NULL;
	xy[1].x = 757;
	xy[1].y = 430;
	xy[1].name = 'A';
	xy[1].exist = true;
	for (i = 2; i <= VerNum; i++)
	{
		xy[i].name = char(i - 1) + 'A';
		xy[i].exist = false;
	}
	for (i = 1; i <= VerNum; i++)
	{
		if (!xy[i].exist)
			continue;
		p = VerList[i].firstEdge;
		while (p)
		{
			if (p->eInfo >= 4 && p->eInfo <= 20)
				generate(xy, i, p);
			p = p->next;
		}
	}
	while (1)
	{
		all = true;
		for (i = 2; i <= VerNum; i++)
		{
			if (!xy[i].exist)
			{
				all = false;
				p = VerList[i].firstEdge;
				while (p)
				{
					if (xy[p->adjVer].exist && p->eInfo >= 4 && p->eInfo <= 20)
						break;
					p = p->next;
				}

				if (p)
				{
					k = p->adjVer;
					p = VerList[k].firstEdge;
					while (p && p->adjVer != i)
						p = p->next;
					generate(xy, k, p);
				}
			}
		}
		if (all)
			break;
	}
	adjustXY(xy);
}

void adjustXY(vN* xy)																//调整结点坐标使分布均匀
{
	int i, j, weight;
	bool all;
	int X = 0, Y = 0;
	while (1)
	{
		all = true;;
		for (i = 1; i < VerNum; i++)
		{
			for (j = i + 1; j <= VerNum; j++)
			{
				getWeight(weight);
				if (sqrt((xy[i].x - xy[j].x) * (xy[i].x - xy[j].x) + (xy[i].y - xy[j].y) * (xy[i].y - xy[j].y)) < 135.0)
				{
					all = false;
					while (sqrt((xy[i].x - xy[j].x) * (xy[i].x - xy[j].x) + (xy[i].y - xy[j].y) * (xy[i].y - xy[j].y)) < 135.0)
					{
						X = rand() % 1375 + 75;
						Y = rand() % 770 + 45;
						xy[j].x = X;
						xy[j].y = Y;
						//cout << "X: " << X << "Y: " << Y << endl;
					}
				}
			}
		}
		if (all)
			break;
	}

}

void paint(vN* xy)																	//显示结点和相关信息
{
	int i;
	EdgeNode* p = NULL;
	setfillcolor(GREEN);
	settextstyle(20, 10, _T("黑体"));
	setbkmode(TRANSPARENT);
	setcolor(BLACK);
	for (i = 1; i <= VerNum; i++)
	{
		fillcircle(xy[i].x, xy[i].y, 8);
		outtextxy(xy[i].x, xy[i].y + 6, xy[i].name);
	}
}

void drawline(vN* xy)																//显示点之间的边
{
	int i;
	EdgeNode* p = NULL;
	setlinecolor(BLACK);
	setlinestyle(PS_SOLID, 1);
	for (i = 1; i <= VerNum; i++)
	{
		p = VerList[i].firstEdge;
		while (p)
		{
			if (p->eInfo >= 4 && p->eInfo <= 20)
				line(xy[i].x, xy[i].y, xy[p->adjVer].x, xy[p->adjVer].y);
			p = p->next;
		}
	}
}

void getd(double& d, double k)														//根据斜率调整步长
{
	k = fabs(k);
	if (k >= 0 && k < 2)
	{
		d = 51.0;
	}
	else if (k >= 2 && k < 5)
	{
		d = 23;
	}
	else if (k >= 5 && k < 10)
	{
		d = 19.0;
	}
	else if (k > 10 && k <= 15)
	{
		d = 5.0;
	}
	else if (k > 15 && k <= 20)
	{
		d = 2.5;
	}
	else if (k > 20 && k <= 30)
	{
		d = 1.75;
	}
	else if (k > 30 && k <= 80)
	{
		d = 1.0;
	}
	else if (k > 80 && k <= 120)
	{
		d = 0.5;
	}
	else if (k > 120 && k <= 250)
	{
		d = 0.35;
	}
	else if (k > 250 && k <= 400)
	{
		d = 0.20;
	}
	else if (k > 400 && k <= 500)
	{
		d = 0.10;
	}
	else if (k > 500 && k <= 800)
	{
		d = 0.08;
	}
	else if (k > 800 && k <= 1200)
	{
		d = 0.05;
	}
	else
	{
		d = 0.01;
	}
}

void flag(int x, int y)																//每到达一点，高亮显示位置
{
	setfillcolor(BROWN);
	fillcircle(x, y, 8);
}

void move(int x, int y, int x1, int y1, vN* xy)										//邮递员移动函数
{
	double d, k, x0, y0, x2, y2;
	x0 = x, y0 = y, x2 = x1, y2 = y1;
	k = (y1 - double(y)) / (x1 - double(x));
	getd(d, k);
	if (x < x1)
	{
		BeginBatchDraw();															//开始批量绘图
		while (x0 < x2)
		{
			paint(xy);
			outtextxy(5, 5, _T("注：图中距离与实际权值并非一一对应，点的相对位置与真实位置未必相符"));
			int j = int(k * (x0 - double(x2)) + y2 + 0.5);
			setfillcolor(DARKGRAY);
			fillcircle(int(x0 + 0.5), j, 2);
			Sleep(500);
			FlushBatchDraw();														//显示
			if (x0 + d >= x2)
			{
				x0 = x2;
				flag(int(x0 + 0.5), int(k * (x0 - double(x2)) + y2 + 0.5));
			}
			else
				x0 += d;
		}
		EndBatchDraw();															   //结束批量绘图
	}
	else
	{
		BeginBatchDraw();
		while (x0 > x2)
		{
			paint(xy);
			outtextxy(5, 5, _T("注：图中距离与实际权值并非一一对应，点的相对位置与真实位置未必相符"));
			setfillcolor(DARKGRAY);
			int j = int(k * (x0 - double(x2)) + y2 + 0.5);
			fillcircle(int(x0 + 0.5), j, 2);
			Sleep(500);
			FlushBatchDraw();
			if (x0 - d <= x2)
			{
				x0 = x2;
				flag(int(x0 + 0.5), int(k * (x0 - double(x2)) + y2 + 0.5));
			}
			else
				x0 -= d;
		}
		EndBatchDraw();
	}
}

void show(vN* xy, Node* Euler)														//动态演示求解过程
{
	Node* p = NULL, * q = NULL;
	p = Euler->next;
	q = p->next;
	while (q)
	{
		move(xy[p->point].x, xy[p->point].y, xy[q->point].x, xy[q->point].y, xy);
		p = q;
		q = q->next;
	}
	paint(xy);
	drawline(xy);
}

int main()
{
	//initgraph(1545, 875, EW_SHOWCONSOLE);											//初始化带有控制台的窗口
	//setbkcolor(WHITE);																//调整背景为白色
	//cleardevice();
	//HWND hwnd = GetHWnd();															//全屏显示
	//SetWindowLong(hwnd, GWL_STYLE, GetWindowLong(hwnd, GWL_STYLE) - WS_CAPTION);
	//SetWindowPos(hwnd, HWND_TOP, 0, 0, GetSystemMetrics(SM_CXSCREEN), GetSystemMetrics(SM_CXSCREEN), SWP_SHOWWINDOW);
	//srand((signed)time(NULL));														//确定随机种子
	initialize();
	vN xy[MaxVerNum];
	Node* Euler = (Node*)malloc(sizeof(Node));
	cout << "Loading a Graph with " << VerNum << " nodes" << endl;
	connected();
	odd();
	Floyd();
	extend_edge();
	search_euler(1);
	get_path(Euler);
	//getXY(xy);
	//show(xy, Euler);
	system("pause");
	//closegraph();																	//关闭窗口
	free();																			//释放内存
	return 0;
}
