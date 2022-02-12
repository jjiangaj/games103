////////////////////////////////////////////////
/// 1. 使用Laplacian smoothing可以明显提高稳定性，smoothing前模型很容易炸。
/// 2. Bonus Task中使用Hyperelastic model会使得模拟速度变慢（在我的电脑上帧率降到大约5，不使用的话大概是20），
/// 得出的模拟效果与不使用Hyperelastic model的StVK没有显著区别。
/// 提交的程序中备注掉了未使用Hyperelastic model的部分。有尝试过Neo-Hookean，不过不太成功，后续再找paper试试看。
/// 3. 调试参数对结果的影响比较大（炸和不炸的区别）。
/// 4. 同学的建议下使用了Parallel.For，提高速度的效果比较显著。
/// ////////////////////////////////////////////////

using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.IO;
using System.Threading.Tasks;

public class FVM : MonoBehaviour
{
	float dt 			= 0.002f;
    float mass 			= 1;
	float stiffness_0	= 7500.0f;
    float stiffness_1 	= 5000.0f;
    float damp			= 0.999f;
    float miu_t = 0.5f;
    float restitution 	= 0.5f;					// for collision


    int[] 		Tet;
	int tet_number;			//The number of tetrahedra

	Vector3[] 	Force;
	Vector3[] 	V;
	Vector3[] 	X;
	Vector3 Fg = new Vector3(0, -9.8f, 0);
	int number;				//The number of vertices

	Matrix4x4[] inv_Dm;

	//For Laplacian smoothing.
	Vector3[]   V_sum;
	int[]		V_num;
	SVD svd = new SVD();

    // Start is called before the first frame update
    void Start()
    {
    	// FILO IO: Read the house model from files.
    	// The model is from Jonathan Schewchuk's Stellar lib.
    	{
    		string fileContent = File.ReadAllText("Assets/house2.ele");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		
    		tet_number=int.Parse(Strings[0]);
        	Tet = new int[tet_number*4];

    		for(int tet=0; tet<tet_number; tet++)
    		{
				Tet[tet*4+0]=int.Parse(Strings[tet*5+4])-1;
				Tet[tet*4+1]=int.Parse(Strings[tet*5+5])-1;
				Tet[tet*4+2]=int.Parse(Strings[tet*5+6])-1;
				Tet[tet*4+3]=int.Parse(Strings[tet*5+7])-1;
			}
    	}
    	{
			string fileContent = File.ReadAllText("Assets/house2.node");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		number = int.Parse(Strings[0]);
    		X = new Vector3[number];
       		for(int i=0; i<number; i++)
       		{
       			X[i].x=float.Parse(Strings[i*5+5])*0.4f;
       			X[i].y=float.Parse(Strings[i*5+6])*0.4f;
       			X[i].z=float.Parse(Strings[i*5+7])*0.4f;
       		}
    		//Centralize the model.
	    	Vector3 center=Vector3.zero;
	    	for(int i=0; i<number; i++)		center+=X[i];
	    	center=center/number;
	    	for(int i=0; i<number; i++)
	    	{
	    		X[i]-=center;
	    		float temp=X[i].y;
	    		X[i].y=X[i].z;
	    		X[i].z=temp;
	    	}
		}
        /*tet_number=1;
        Tet = new int[tet_number*4];
        Tet[0]=0;
        Tet[1]=1;
        Tet[2]=2;
        Tet[3]=3;

        number=4;
        X = new Vector3[number];
        V = new Vector3[number];
        Force = new Vector3[number];
        X[0]= new Vector3(0, 0, 0);
        X[1]= new Vector3(1, 0, 0);
        X[2]= new Vector3(0, 1, 0);
        X[3]= new Vector3(0, 0, 1);*/


        //Create triangle mesh.
       	Vector3[] vertices = new Vector3[tet_number*12];
        int vertex_number=0;
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];

        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }

        int[] triangles = new int[tet_number*12];
        for(int t=0; t<tet_number*4; t++)
        {
        	triangles[t*3+0]=t*3+0;
        	triangles[t*3+1]=t*3+1;
        	triangles[t*3+2]=t*3+2;
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.triangles = triangles;
		mesh.RecalculateNormals ();


		V 	  = new Vector3[number];
        Force = new Vector3[number];
        V_sum = new Vector3[number];
        V_num = new int[number];

		//TODO: Need to allocate and assign inv_Dm
		inv_Dm = new Matrix4x4[tet_number];
		for(int t=0; t<tet_number; t++)
		{
			inv_Dm[t] = Build_Edge_Matrix(t).inverse;
		}

		for (int i = 0; i < number; i++)
		{
			V_sum[i] = new Vector3(0, 0, 0);
			V_num[i] = 0;
		}
    }

    Matrix4x4 Build_Edge_Matrix(int tet)
    {
    	Matrix4x4 ret=Matrix4x4.zero;
    	//TODO: Need to build edge matrix here.
        Vector3 X0 = X[Tet[tet*4+0]];
        Vector3 X10 = X[Tet[tet*4+1]] - X0;
        Vector3 X20 = X[Tet[tet*4+2]] - X0;
        Vector3 X30 = X[Tet[tet*4+3]] - X0;
        ret.SetColumn(0, X10);
        ret.SetColumn(1, X20);
        ret.SetColumn(2, X30);
        ret[3, 3] = 1;
        return ret;
    }
    
    Matrix4x4 MatAdd(Matrix4x4 m, Matrix4x4 n)
    {
	    Matrix4x4 A = Matrix4x4.zero;
	    for(int i = 0; i < 4; i++){
		    for(int j = 0; j < 4; j++){
			    A [i, j] =  m [i,j] + n[i, j];
		    }
	    }
	    A[3, 3] = 1;
	    return A;
    }

    Matrix4x4 MatSub(Matrix4x4 m, Matrix4x4 n)
    {
	    Matrix4x4 A = Matrix4x4.zero;
	    for(int i = 0; i < 4; i++){
		    for(int j = 0; j < 4; j++){
			    A [i, j] =  m [i,j] - n[i, j];
		    }
	    }
	    A[3, 3] = 1;
	    return A;
    }
    
    Matrix4x4 Times(Matrix4x4 m, float k)
    {
	    Matrix4x4 M = m;
	    for(int i = 0; i < 4; i++){
		    for(int j = 0; j < 4; j++){
			    M [i, j] *= k;
		    }
	    }
	    M[3, 3] = 1;
	    return M;
    }


    void _Update()
    {
    	// Jump up.
		if(Input.GetKeyDown(KeyCode.Space))
    	{
    		for(int i=0; i<number; i++)
    			V[i].y+=0.3f;
    	}

    	for(int i=0 ;i<number; i++)
    	{
    		//TODO: Add gravity to Force.
	        Force[i] = Fg;
        }
        
        //Without Hyper elastic-----------------------------------------///
     //    Parallel.For(0, tet_number, (tet) =>
    	// // for(int tet=0; tet<tet_number; tet++)
    	// {
    	// 	//TODO: Deformation Gradient
     //        Matrix4x4 dm = Build_Edge_Matrix(tet);
     //        Matrix4x4 F = dm * inv_Dm[tet];
     //  
     //        //TODO: Green Strain
     //        Matrix4x4 G = Times(MatSub(F.transpose * F, Matrix4x4.identity), 0.5f);
     //        float trace = G[0, 0] + G[1, 1] + G[2, 2];
     //
     //        //TODO: First PK Stress
     //        Matrix4x4 P = F * MatAdd(Times(G, 2 * stiffness_1), Times(Matrix4x4.identity, stiffness_0 * trace));
     //  
     //        //TODO: Elastic Force
     //        float param = -1 / (6 * inv_Dm[tet].determinant);
     //        Matrix4x4 Fi = Times(P, param) * inv_Dm[tet].transpose;
     //  
     //        Force[Tet[tet * 4 + 0]] -= (Vector3) (Fi.GetColumn(0) + Fi.GetColumn(1) + Fi.GetColumn(2));
     //        Force[Tet[tet * 4 + 1]] += (Vector3) Fi.GetColumn(0);
     //        Force[Tet[tet * 4 + 2]] += (Vector3) Fi.GetColumn(1);
     //        Force[Tet[tet * 4 + 3]] += (Vector3) Fi.GetColumn(2);
     //    });
		//--------------------------------------------------------------///
		
		
		//Hyperelastic--------------------------------------------------///
        Parallel.For(0, tet_number, (tet) =>
	        // for(int tet=0; tet<tet_number; tet++)
        {
	        //TODO: Deformation Gradient
	        Matrix4x4 dm = Build_Edge_Matrix(tet);
	        Matrix4x4 F = dm * inv_Dm[tet];
	        Matrix4x4 C = F.transpose*F;

	        //SVD
	        float I = C[0, 0] + C[1, 1] + C[2, 2];
	        
	        Matrix4x4 U = Matrix4x4.zero;
	        Matrix4x4 S = Matrix4x4.zero;
	        Matrix4x4 V = Matrix4x4.zero;
	        svd.svd(F, ref U, ref S, ref V);
	        
	        Matrix4x4 diag = Matrix4x4.zero;
	        diag[3, 3] = 1;
	        diag[0, 0] = stiffness_0 * (I - 3) * S[0,0]/2 + stiffness_1 * (S[0,0]*S[0,0]*S[0,0] - S[0,0]);
	        diag[1, 1] = stiffness_0 * (I - 3) * S[1,1]/2 + stiffness_1 * (S[1,1]*S[1,1]*S[1,1] - S[1,1]);
	        diag[2, 2] = stiffness_0 * (I - 3) * S[2,2]/2 + stiffness_1 * (S[2,2]*S[2,2]*S[2,2] - S[2,2]);
	        Matrix4x4 P = U * diag * V.transpose;
      
	        //TODO: Elastic Force
	        float param = -1 / (6 * inv_Dm[tet].determinant);
	        Matrix4x4 Fi = Times(P, param) * inv_Dm[tet].transpose;
      
	        Force[Tet[tet * 4 + 0]] -= (Vector3) (Fi.GetColumn(0) + Fi.GetColumn(1) + Fi.GetColumn(2));
	        Force[Tet[tet * 4 + 1]] += (Vector3) Fi.GetColumn(0);
	        Force[Tet[tet * 4 + 2]] += (Vector3) Fi.GetColumn(1);
	        Force[Tet[tet * 4 + 3]] += (Vector3) Fi.GetColumn(2);
        });
        //--------------------------------------------------------------///
        
        //Laplacian smoothing
        Parallel.For(0, tet_number, (tet) =>
	        // for(int tet=0; tet<tet_number; tet++)
        {
	        //For a tetrahedral, add the 3 velocities to the top vertex
	        Vector3 v0 = V[Tet[tet * 4 + 0]];
	        Vector3 v1 = V[Tet[tet * 4 + 1]];
	        Vector3 v2 = V[Tet[tet * 4 + 2]];
	        Vector3 v3 = V[Tet[tet * 4 + 3]];
	        Vector3 sum = v0 + v1 + v2 + v3;

	        //Account for the 3 edges linked to the top vertex
	        V_sum[Tet[tet * 4 + 0]] = sum - v0;
	        V_num[Tet[tet * 4 + 0]] = 3;

	        //Add the top vertex's velocity to the other vertices
	        V_sum[Tet[tet * 4 + 1]] = sum - v1;
	        V_num[Tet[tet * 4 + 1]] = 3;

	        V_sum[Tet[tet * 4 + 2]] = sum - v2;
	        V_num[Tet[tet * 4 + 2]] = 3;

	        V_sum[Tet[tet * 4 + 3]] = sum - v3;
	        V_num[Tet[tet * 4 + 3]] = 3;
        });

        Parallel.For(0, number, (i) =>
	        // for (int i = 0; i < number; i++)
        {
	        V[i] = 0.95f * V[i] + 0.05f * (V_sum[i] / V_num[i]);
        });

        for (int i = 0; i < number; i++)
        {
	        //TODO: Update V and X here。
	        V[i] += dt * Force[i];
	        V[i] *= damp;
	        X[i] += dt * V[i];
        }

        //TODO: (Particle) collision with floor.
        Vector3 P1 = new Vector3(0, -3, 0);
        Vector3 N1 = new Vector3(0, 1, 0);
        Parallel.For(0, number, (i) =>
	        // for (int i = 0; i < number; i++)
        {
	        float phi = Vector3.Dot((X[i] - P1), N1);

	        if (phi < 0)
	        {
		        X[i] = X[i] - phi * N1;
		        if (Vector3.Dot(V[i], N1) < 0)
		        {
			        Vector3 vn = Vector3.Dot(V[i], N1) * N1;
			        Vector3 vt = V[i] - vn;

			        float a = Mathf.Max(1 - miu_t * (1 + restitution) * vn.magnitude / vt.magnitude, 0);
			        // 1 - μt(1 + μn)||Vni||/||Vti||

			        vn *= -restitution;
			        vt *= a;
			        V[i] = vn + vt;
		        }
	        }
        });
        
        
    }
	
    

    // Update is called once per frame
    void Update()
    {
    	for(int l=0; l<10; l++)
    		 _Update();

    	// Dump the vertex array for rendering.
    	Vector3[] vertices = new Vector3[tet_number*12];
        int vertex_number=0;
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.RecalculateNormals ();
    }
}
