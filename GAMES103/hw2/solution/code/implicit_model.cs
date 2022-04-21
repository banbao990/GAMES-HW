using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class implicit_model : MonoBehaviour {
    const float t = 0.0333f;
    const float t_neg = 1 / t;
    const float t2_neg = 1 / (t * t);
    const float mass = 1;
    const float damping = 0.99f;
    const float rho = 0.8F; // 0.995f 切比雪夫会炸
    const float spring_k = 8000;
    const float spring_k4 = spring_k * 4;
    static readonly HashSet<int> fixedPoint = new HashSet<int> { 0, 20 };
    const int N = 21;       // 将 mesh 重构为 20*20 的网格

    int[] E;                // 弹簧边的信息(结点索引), 都是直的弹簧, 不能被弯曲
    float[] L;              // 每条边的原长
    Vector3[] V;            // V 结点的速度信息
                            // X 是结点位置信息, 节点 0/20 为固定点(更新的时候直接忽略即可)
                            // 结点数: 21 * 21 = 441

    bool useChebyshev = false;

    // Start is called before the first frame update
    void Start() {
        GameObject.Find("Canvas").GetComponent<Canvas>().enabled = true;

        Mesh mesh = GetComponent<MeshFilter>().mesh;

        // Resize the mesh.

        // 更新到 mesh 的变量
        Vector3[] X = new Vector3[N * N];
        Vector2[] UV = new Vector2[N * N]; // 纹理坐标
        int[] triangles = new int[(N - 1) * (N - 1) * 6]; // 三个一组表示一个三角形, 数量 20*20*2(*3 结点)
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                X[j * N + i] = new Vector3(5 - 10.0f * i / (N - 1), 0, 5 - 10.0f * j / (N - 1)); // 范围 [5, -5]*[5, -5], i 对应 z
                                                                                                 // x: 右, y: 上, z: 里
                UV[j * N + i] = new Vector2(i / (N - 1.0f), j / (N - 1.0f));
            }
        }
        int t = 0;
        for (int j = 0; j < N - 1; j++) {
            for (int i = 0; i < N - 1; i++) {
                triangles[t * 6 + 0] = j * N + i;
                triangles[t * 6 + 1] = j * N + i + 1;
                triangles[t * 6 + 2] = (j + 1) * N + i + 1;
                triangles[t * 6 + 3] = j * N + i;
                triangles[t * 6 + 4] = (j + 1) * N + i + 1;
                triangles[t * 6 + 5] = (j + 1) * N + i;
                t++;
            }
        }
        mesh.vertices = X;
        mesh.triangles = triangles;
        mesh.uv = UV;
        mesh.RecalculateNormals();

        // Construct the original E
        int[] _E = new int[triangles.Length * 2];
        // 每个三角形拆分为 3 条边(相邻的两个结点 [2i, 2i+1] 表示一条边)
        for (int i = 0; i < triangles.Length; i += 3) {
            _E[i * 2 + 0] = triangles[i + 0];
            _E[i * 2 + 1] = triangles[i + 1];
            _E[i * 2 + 2] = triangles[i + 1];
            _E[i * 2 + 3] = triangles[i + 2];
            _E[i * 2 + 4] = triangles[i + 2];
            _E[i * 2 + 5] = triangles[i + 0];
        }

        // Reorder the original edge list
        // 使得第一个顶点序号比第二个小(排序之后去重)
        for (int i = 0; i < _E.Length; i += 2)
            if (_E[i] > _E[i + 1])
                Swap(ref _E[i], ref _E[i + 1]);

        // Sort the original edge list using quicksort
        Quick_Sort(ref _E, 0, _E.Length / 2 - 1);

        // 去重(统计边数)
        int e_number = 0;
        for (int i = 0; i < _E.Length; i += 2) {
            if (i == 0 || _E[i + 0] != _E[i - 2] || _E[i + 1] != _E[i - 1]) {
                e_number++;
            }
        }

        // 去重, 拷贝到 E
        E = new int[e_number * 2];
        for (int i = 0, e = 0; i < _E.Length; i += 2) {
            if (i == 0 || _E[i + 0] != _E[i - 2] || _E[i + 1] != _E[i - 1]) {
                E[e * 2 + 0] = _E[i + 0];
                E[e * 2 + 1] = _E[i + 1];
                e++;
            }
        }

        // 弹簧原长
        L = new float[E.Length / 2];
        for (int e = 0; e < E.Length / 2; e++) {
            int v0 = E[e * 2 + 0];
            int v1 = E[e * 2 + 1];
            L[e] = (X[v0] - X[v1]).magnitude;
        }

        V = new Vector3[X.Length];
        for (int i = 0; i < V.Length; i++) {
            V[i] = Vector3.zero;
        }
    }

    void Collision_Handling() {
        // [1.e] Sphere Collision
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] X = mesh.vertices;
        GameObject sphere = GameObject.Find("Sphere");
        Vector3 center = sphere.transform.position;
        const float r = 2.7F;
        const float r2 = r * r;

        // Handle collision
        for (int i = 0; i < X.Length; ++i) {
            if (fixedPoint.Contains(i)) { continue; }
            Vector3 disV = X[i] - center;
            float dis2 = disV.sqrMagnitude;
            if (dis2 < r2) {
                // collision
                float dis = Mathf.Sqrt(dis2);
                Vector3 disR = r * disV / dis;
                V[i] += t_neg * (center + disR - X[i]);
                X[i] = center + disR;
            }
        }

        mesh.vertices = X;
    }

    void Get_Gradient(Vector3[] X, Vector3[] X_hat, float t, Vector3[] G) {
        // [1.b] Gradient Calculation

        int length = G.Length;
        Vector3 gravity_neg = mass * new Vector3(0, 9.8F, 0);

        // 注意正负, 这里是梯度
        // Momentum and Gravity.
        for (int i = 0; i < length; ++i) {
            G[i] = t2_neg * mass * (X[i] - X_hat[i]) + gravity_neg;
        }

        // Spring Force.
        length = E.Length / 2;
        for (int i = 0; i < length; ++i) {
            int tmp = 2 * i;
            int l = E[tmp], r = E[tmp + 1];
            Vector3 dis = X[l] - X[r];
            Vector3 f_neg = spring_k * (1 - L[i] / dis.magnitude) * dis;

            G[l] += f_neg;
            G[r] -= f_neg;
        }

        foreach(int i in fixedPoint) {
            G[i] = Vector3.zero;
        }
    }

    // Update is called once per frame
    void Update() {
        // Game Control
        if (Input.GetKey("s")) {
            ResetCloth();
        }
        if (Input.GetKeyDown(KeyCode.Space)) {
            useChebyshev = !useChebyshev;
            Text info = GameObject.Find("MethodInfo").GetComponent<Text>();
            if (useChebyshev) {
                info.text = "切比雪夫方法";
            } else {
                info.text = "类牛顿法";
            }
        }

        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] X = mesh.vertices;
        Vector3[] last_X = new Vector3[X.Length];
        Vector3[] X_hat = new Vector3[X.Length];
        Vector3[] G = new Vector3[X.Length];

        int N = X.Length;

        // [1.a] Initial Setup.
        for (int i = 0; i < N; ++i) {
            if (fixedPoint.Contains(i)) {
                V[i] = Vector3.zero;
                X_hat[i] = X[i];
                continue;
            }

            V[i] *= damping;
            X_hat[i] = X[i] + t * V[i];
        }

        // initial guess: 猜测的好坏会影响收敛速度
        for (int i = 0; i < N; ++i) {
            X[i] = X_hat[i];
        }


        if (useChebyshev) {
            // [1.d] Chebyshev Acceleration
            float omega = 1.0F;
            const float rho2 = rho * rho;
            float A_neg = 1 / (t2_neg * mass + spring_k4);

            for (int k = 0; k < 32; k++) {
                if (k == 0) {
                    omega = 1.0F;
                } else if (k == 1) {
                    omega = 2 / (2 - rho2);
                } else {
                    omega = 4 / (4 - rho2 * omega);
                }
                float one_minus_omega = 1 - omega;

                Get_Gradient(X, X_hat, t, G);

                // 切比雪夫加速是因为原来求 A^-1 很麻烦才提出的, 但是 A^-1 很容易求 
                // 这里只是用了 omega*deltaX + (1-omega)*last_deltaX 的思想
                // 等价于       omega*X      + (1-omega)*last_X

                for (int i = 0; i < N; ++i) {
                    if (fixedPoint.Contains(i)) { continue; }
                    Vector3 new_X = omega * (X[i] - A_neg * G[i]) + one_minus_omega * last_X[i];
                    last_X[i] = X[i];
                    X[i] = new_X;
                }
            }
        } else {
            // [1.c] Update X by gradient.
            float A_neg = 1 / (t2_neg * mass + spring_k4);
            for (int k = 0; k < 32; k++) {
                Get_Gradient(X, X_hat, t, G);
                for (int i = 0; i < N; ++i) {
                    if (fixedPoint.Contains(i)) { continue; }
                    X[i] -= A_neg * G[i];
                }
            }
        }

        // [1.c] Finishing
        // 更新位置和速度
        for (int i = 0; i < N; ++i) {
            V[i] += t_neg * (X[i] - X_hat[i]);
        }

        mesh.vertices = X;

        Collision_Handling();
        mesh.RecalculateNormals();
    }

    void ResetCloth() {
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] X = mesh.vertices;
        for (int i = 0; i < X.Length; ++i) {
            V[i] = Vector3.zero;
        }
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                X[j * N + i] = new Vector3(5 - 10.0f * i / (N - 1), 0, 5 - 10.0f * j / (N - 1));
            }
        }
        mesh.vertices = X;
    }

    #region QuickSort Module
    void Quick_Sort(ref int[] a, int l, int r) {
        int j;
        if (l < r) {
            j = Quick_Sort_Partition(ref a, l, r);
            Quick_Sort(ref a, l, j - 1);
            Quick_Sort(ref a, j + 1, r);
        }
    }

    int Quick_Sort_Partition(ref int[] a, int l, int r) {
        // a[2*i], a[2*i+1] 是一个整体
        int pivot_0, pivot_1, i, j;
        pivot_0 = a[l * 2 + 0];
        pivot_1 = a[l * 2 + 1];
        i = l;
        j = r + 1;
        while (true) {
            do ++i; while (i <= r && (a[i * 2] < pivot_0 || a[i * 2] == pivot_0 && a[i * 2 + 1] <= pivot_1));
            do --j; while (a[j * 2] > pivot_0 || a[j * 2] == pivot_0 && a[j * 2 + 1] > pivot_1);
            if (i >= j) break;
            Swap(ref a[i * 2], ref a[j * 2]);
            Swap(ref a[i * 2 + 1], ref a[j * 2 + 1]);
        }
        Swap(ref a[l * 2 + 0], ref a[j * 2 + 0]);
        Swap(ref a[l * 2 + 1], ref a[j * 2 + 1]);
        return j;
    }

    void Swap(ref int a, ref int b) {
        int temp = a;
        a = b;
        b = temp;
    }

    #endregion QuickSort Module
}
