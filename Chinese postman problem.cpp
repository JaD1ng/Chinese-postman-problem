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

/*���˼�룺
���ѡȡ����ߣ�ͬʱ�������鼯����ͼ����ͨ�ԣ������ɵ�ͼ����ͨ����ߣ����Բ����������������¹���ֱ������������
ͨ��Floyd�㷨��þ��������·�����飬���ں���ƥ��Ͳ��ҡ��������±��Լ�����������״ѹdp������СȨƥ�䣬������ԡ�
��ԭͼ����ӱ����������ŷ��ͼ����DFS�ҵ�һ��ŷ����·����ŷ����·�У����Ʋ�ֱ����������·�������ʵ��·����
��ȡ������꣬�̶��ʾ�A�㣬���ݽ�����ȷ��Ȩֵ���ж��������룬���̫С������ѡȡ�����յõ��ֲ����ȵ�ͼ��
�������н�㣬��㰴·���ƶ����������ߣ���̬��ʾ��ϡ�
*/

#define INF 65535			//���������
#define MaxVerNum  100		//�������Ŀ
#define eps 1e-3			//���ȣ��ھ���֮�������Ҫ��

string num_to_uppercase(int num) {
    // �������ֺ���ĸ�Ķ�Ӧ��ϵ
    char num_to_letter[] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};
    // ��������Ƿ��� 1-26 ��Χ�ڣ����򷵻ؿ��ַ���
    if (num < 1 || num > 26) {
        return "";
    }
    // ������ת��Ϊ��д��ĸ
    return string(1, num_to_letter[num - 1]);
}

typedef char elementType;	//����ͼ�ж������������
typedef int eInfoType;		//���� eInfo ���������ͣ���Ȩֵ����������

typedef struct eNode		//���������Ľ��ṹ
{
	int adjVer;				//�ڽӶ�����Ϣ���˴�Ϊ�����ţ��� 1 ��ʼ
	eInfoType eInfo;		//�������б�ʾ�ߵ������Ϣ��������Ȩֵ
	struct eNode* next;		//ָ��������е���һ����㡣
	bool visited;
}EdgeNode;					//������������

typedef struct vNode		//���嶥���Ľ��ṹ
{
	elementType data;		//�������������
	EdgeNode* firstEdge;	//ָ��߽��
}VerNode;					//�����������

typedef struct oNode {		//���Դ洢�ṹ
	int node_1;
	int node_2;
}OddNode;

typedef struct Euler {		//·����㶨��
	int point;				//����±�
	struct Euler* next;		//ָ����һ��·�����
}Node;

typedef  struct {			//�����������Ϣ�Ľṹ
	int x, y;				//�ᡢ������
	elementType name;		//�������A ~ Z
	bool exist;				//�Ƿ��ѻ������
}vN;

stack<int> Q;													//ȫ��ջ�����ڻ�ý��֮���·��
int bg, to;														//�ߵ������յ�
int VerNum;														//������
int odd_num = 0, lf = 0;										//������
int* dp;														//ȫ��ָ�룬������������̬�����ڴ�
bool isIP = true;												//�����ж�VerNum�Ƿ�����ֵ
int degree[MaxVerNum];											//��ʾ������
int tree[MaxVerNum];											//���鼯�������ж�ͼ����ͨ��
int path[MaxVerNum][MaxVerNum];									//·������
VerNode VerList[MaxVerNum];										//�����
eInfoType dist[MaxVerNum][MaxVerNum];							//�������
eInfoType AdjMatrix[MaxVerNum][MaxVerNum];						//�ڽӾ���
int* odd_point;													//��¼���
int odd_n[MaxVerNum];											//���ڼ�¼������������ǵڼ������
OddNode odd_match[MaxVerNum];									//�������
void initialize();												//��ʼ������
int findroot(int x);											//�������鼯
bool add(int m, int n, int w, bool instant);					//��ӱ�
void free();													//�ͷŶ����ı�������̬����������
void connected();												//�����ɵ�ͼ����ͨ����ϲ�
bool Legal();													//�ж�ͼ�Ƿ��������
void odd();														//��������������������ͼ
void Floyd();													//������������������Լ�·��
int get_oddedge(int x);											//�����Ҫ��ߵ�����
void extend_edge();												//��ԭͼ����ӱ��Ի��ŷ��ͼ
bool IsExtend(int x, int y);									//�ж��Ƿ�����ӵı�
void search_euler(int u);										//Ѱ��һ��ŷ����·����¼
void get_path(Node* p);											//���ʵ��·��
void getarc(int& arc, int x, int y);							//���ݵ��λ�û�ú��ʵĽǶ�
void getWeight(int& weight);									//���ݵ����Ŀȷ�������ľ���
void getXY(vN* xy);												//��ȡ�������
void generate(vN* xy, int i, EdgeNode* p);						//�����������
void adjustXY(vN* xy);											//����������꣬ʹ����֮�䱣��һ������
void paint(vN* xy);												//���Ƴ����н��
void drawline(vN* xy);											//���Ƶ�֮��ı�
void show(vN* xy, Node* Euler);									//��̬��ʾ�����̣�·����
void getd(double& d, double k);									//��б�ʲ�ͬѡȡ���ʵĲ���
void flag(int x, int y);										//����һ����������һ�Σ�����ȷ���ʵ�Ա��λ��
void move(int x1, int y1, int x2, int y2, vN* xy);				//��ʾ�ʵ�Ա��ͼ�ϵ��ƶ�����

void initialize()																	//��ʼ������
{
	int count = 0, i, m, n, weight, EdgeNum;
	bool succeed;
	if (isIP)																		//�Ƿ���Ҫ�ٴ����ɽ�����
	{
		VerNum = rand() % 16 + 9;													//��������9 ~ 24
		isIP = false;
	}
	for (i = 1; i < VerNum + 1; i++)												//��ʼ�������
	{
		VerList[i].data = char(i - 1) + 'A';
		VerList[i].firstEdge = NULL;
	}
	odd_point = (int*)malloc((VerNum + 1) * sizeof(int));						//��̬����
	memset(degree, 0, sizeof(degree));									//��ʼ��dgree����
	memset(tree, -1, sizeof(tree));									//��ʼ��tree����
	EdgeNum = rand() % (VerNum * 2 - VerNum + 1) + VerNum - 1;						//���Ʊߵĸ���
	while (count < EdgeNum)															//ȷ��������Լ�Ȩֵ
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

int findroot(int x)																	//���鼯�ķ�ʽ�ж�ͼ��ͨ��
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

bool add(int m, int n, int weight, bool instant)									//��ߣ�����instant��ֵ��ȷ��������ʽ�����
{
	EdgeNode* p = NULL, * edge1 = NULL, * edge2 = NULL;
	int t, tempa, tempb;
	if (!instant)																	//�������ظ��ߵ���߷�ʽ
	{																				//�����������򷵻�false			
		if (degree[m] > 3 && degree[n] > 3 && m < 1 || m > VerNum || n < 1 || n > VerNum)
			return false;
		p = VerList[m].firstEdge;													//����
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
		edge1 = new EdgeNode;														//���˫���
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
		tempa = findroot(m);														//�������鼯�����ж���ͨ��
		tempb = findroot(n);
		if (tempa != tempb)															//�����Ȳ�ͬ��ϲ���ͼ
			tree[tempa] = tempb;
	}
	else                                                                            //�ظ��ߵ����
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
	degree[m]++;																	//��������Ӧ����
	degree[n]++;
	return true;
}

void connected()
{
	int i, ans = 0, a[80], weight;
	memset(a, -1, sizeof(a));											//��ʼ�����Դ�Ŷ�����ͼ���Ƚ�������
	for (i = 1; i < VerNum + 1; i++)
	{																				//�ж���ͨ��
		if (tree[i] == -1)
		{
			a[ans] = i;
			ans++;
		}
	}
	i = 0;
	if (ans > 1)																	//�ϲ���ͼ��������һ����ͨͼ
	{
		while (i < ans - 1)
		{
			weight = rand() % 15 + 4;
			add(a[i], a[i + 1], weight, false);
			i++;
		}
	}
}

bool Legal()																		//�ж���ͼ�Ƿ�Ϸ�
{
	int i, odd = 0;
	for (i = 1; i < VerNum + 1; i++)												//ͳ��������
	{
		if (degree[i] % 2)
			odd++;
	}

	if (odd == 0 || odd > 2)
		return true;
	if (odd == 2)
	{
		if (degree[1] % 2 == 1)														//��������Ϊ2����1��Ϊ��㣬�򲻷���Ҫ��
			return true;
		else
			return false;
	}
}

void odd()																			//ѭ������ͼ��ֱ�����ɵ�ͼ����Ҫ��
{
	while (!Legal())
	{
		free();
		initialize();
		connected();
	}
}

void Floyd()																		//Floyd�㷨ȷ������������뼰·��
{
	int i, j, k;
	EdgeNode* p = NULL;
	for (i = 1; i < VerNum + 1; i++)
		for (j = 1; j < VerNum + 1; j++)
			AdjMatrix[i][j] = INF;
	for (i = 1; i < VerNum + 1; i++)												//���ڽӱ�ת��Ϊ�ڽӾ���
	{
		p = VerList[i].firstEdge;
		while (p)
		{
			AdjMatrix[i][p->adjVer] = p->eInfo;
			p = p->next;
		}
	}
	for (i = 1; i < VerNum + 1; i++)												//��ʼ��·������
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

int get_oddedge(int x)																//�������СȨƥ�������
{
	int i, j;
	int ans = INF, temp, ANS;
	EdgeNode* point = NULL;
	if (dp[x] || !x)
		return dp[x];
	for (i = 0; i < odd_num - 1; i++)												//����ѹ��״̬�е�����'1'�ԣ������е�����
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
						ans = ANS;													//��¼Ȩֵ��С�������±�
						bg = odd_n[odd_num - i];
						to = odd_n[odd_num - j];
					}
				}
			}
		}
	}
	return dp[x] = ans;
}

void extend_edge()																	//������Եı�
{
	int i, j = 0, a = 0, count = 1, n = 0;
	bool instant = false;
	EdgeNode* point = NULL;
	cout << "����У�" << endl;
	for (i = 1; i < VerNum + 1; i++)												//��ʼ���������
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
	cout << "������Ϊ��" << lf << endl;
	for (i = odd_num; i >= 1; i--)													//�������������ɶ�Ӧ��������
	{
		a <<= 1;
		a++;
	}
	i = 0;
	while (a)
	{
		dp = (int*)malloc((1 << odd_num) * sizeof(int));							//��������ʼ��dp����
		memset(dp, 0, (1 << odd_num) * sizeof(int));
		instant = false;
		get_oddedge(a);																//��ȡ���Բ���¼
		odd_match[i].node_1 = bg;
		odd_match[i].node_2 = to;
		i++;
		for (j = 1; j <= odd_num; j++)												//����������飨��ߺ�����ɾ����
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
		point = VerList[bg].firstEdge;												//���
		while (point)																//�ж�������Ƿ�ֱ������
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
		free(dp);																	//�ͷ�dp����
	}
	for (int i = 1; i <= VerNum; i++)
	{
		if (degree[i] % 2)
			n++;
	}
	cout << "��ߺ��������Ϊ��" << n << endl;
}

bool IsExtend(int x, int y)															//�Ƿ�Ϊ������ӵı�
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

void search_euler(int u)															 //DFSѰ��ŷ����·
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
			while (q && q->visited)													//��ֹ�ر��޷���ȷ����
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
	Q.push(u);																		//�ڽӵ㶼����������ջ
}

void get_path(Node* p)																//��ÿ���·��
{
	Node* q = p, * t = NULL, * r = NULL;
	int x, k;
	bool instant = true;
	while (!Q.empty())																//����ŷ����·
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
		cout << "�����޽�" << endl;
		exit(1);
	}
	cout << "ŷ����·Ϊ��" << endl;
	while (r)
	{
		cout << char(r->point - 1 + 'A');
		if (r->next)
			cout << "->";
		r = r->next;
	}
	cout << endl;
	while (t)																	//�õ�����·��
	{
		instant = true;
		if (IsExtend(q->point, t->point))										//��ȫ��ֱ�����������·��
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
	cout << "����·��Ϊ��" << endl;
	while (q)
	{
		cout << char(q->point - 1 + 'A');
		if (q->next)
			cout << "->";
		q = q->next;
	}
	cout << endl;
}

void free()																			//�ͷŶ�̬������ڴ�
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

void getarc(int& arc, int x, int y)													//�������ýǶȣ�ʹ�ֲ�����
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

void getWeight(int& weight)															 //���ݵ�ĸ���ȷ������
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

void generate(vN* xy, int i, EdgeNode* p)											//����ĳ��ΪԲ�ĵ�Բ��ȡ��һ��
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

void getXY(vN* xy)																	//��ý��XY����
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

void adjustXY(vN* xy)																//�����������ʹ�ֲ�����
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

void paint(vN* xy)																	//��ʾ���������Ϣ
{
	int i;
	EdgeNode* p = NULL;
	setfillcolor(GREEN);
	settextstyle(20, 10, _T("����"));
	setbkmode(TRANSPARENT);
	setcolor(BLACK);
	for (i = 1; i <= VerNum; i++)
	{
		fillcircle(xy[i].x, xy[i].y, 8);
		outtextxy(xy[i].x, xy[i].y + 6, xy[i].name);
	}
}

void drawline(vN* xy)																//��ʾ��֮��ı�
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

void getd(double& d, double k)														//����б�ʵ�������
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

void flag(int x, int y)																//ÿ����һ�㣬������ʾλ��
{
	setfillcolor(BROWN);
	fillcircle(x, y, 8);
}

void move(int x, int y, int x1, int y1, vN* xy)										//�ʵ�Ա�ƶ�����
{
	double d, k, x0, y0, x2, y2;
	x0 = x, y0 = y, x2 = x1, y2 = y1;
	k = (y1 - double(y)) / (x1 - double(x));
	getd(d, k);
	if (x < x1)
	{
		BeginBatchDraw();															//��ʼ������ͼ
		while (x0 < x2)
		{
			paint(xy);
			outtextxy(5, 5, _T("ע��ͼ�о�����ʵ��Ȩֵ����һһ��Ӧ��������λ������ʵλ��δ�����"));
			int j = int(k * (x0 - double(x2)) + y2 + 0.5);
			setfillcolor(DARKGRAY);
			fillcircle(int(x0 + 0.5), j, 2);
			Sleep(500);
			FlushBatchDraw();														//��ʾ
			if (x0 + d >= x2)
			{
				x0 = x2;
				flag(int(x0 + 0.5), int(k * (x0 - double(x2)) + y2 + 0.5));
			}
			else
				x0 += d;
		}
		EndBatchDraw();															   //����������ͼ
	}
	else
	{
		BeginBatchDraw();
		while (x0 > x2)
		{
			paint(xy);
			outtextxy(5, 5, _T("ע��ͼ�о�����ʵ��Ȩֵ����һһ��Ӧ��������λ������ʵλ��δ�����"));
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

void show(vN* xy, Node* Euler)														//��̬��ʾ������
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
	//initgraph(1545, 875, EW_SHOWCONSOLE);											//��ʼ�����п���̨�Ĵ���
	//setbkcolor(WHITE);																//��������Ϊ��ɫ
	//cleardevice();
	//HWND hwnd = GetHWnd();															//ȫ����ʾ
	//SetWindowLong(hwnd, GWL_STYLE, GetWindowLong(hwnd, GWL_STYLE) - WS_CAPTION);
	//SetWindowPos(hwnd, HWND_TOP, 0, 0, GetSystemMetrics(SM_CXSCREEN), GetSystemMetrics(SM_CXSCREEN), SWP_SHOWWINDOW);
	//srand((signed)time(NULL));														//ȷ���������
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
	//closegraph();																	//�رմ���
	free();																			//�ͷ��ڴ�
	return 0;
}
