var js3d = {
    vector2 : {
        normalize(v) {
            var l = v[0]*v[0]+v[1]*v[1]
            if (l==0) return v
            l = Math.sqrt(l)
            return [
                v[0] / l,
                v[1] / l,
            ]
        },
        dis(v) {
            var l = v[0]*v[0]+v[1]*v[1]
            if (l==0) return 0
            return Math.sqrt(l)
        },
        clone(v) {
            return [
                v[0],
                v[1],
            ]
        },
        sub(v1,v2) {
            return [
                v1[0] - v2[0],
                v1[1] - v2[1],
            ]
        },
        dot(v1,v2) {
            return v1[0]*v2[0]+v1[1]*v2[1]
        },
        cross(v1, v2) {
            var a=v1[0]
            var b=v1[1]
            var d=v2[0]
            var e=v2[1]
            return a*e-b*d
        },
    },
    vector3 : {
        normalize(v) {
            var l = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]
            if (l==0) return v
            l = Math.sqrt(l)
            return [
                v[0] / l,
                v[1] / l,
                v[2] / l,
            ]
        },
        clone(v) {
            return [
                v[0],
                v[1],
                v[2],
            ]
        },
        disAB(v1,v2) {
            return this.dis(this.sub(v1,v2))
        },
        dis(v) {
            var l = v[0]*v[0]+v[1]*v[1]+v[2]*v[2]
            if (l==0) return 0
            return Math.sqrt(l)
        },
        dot(v1,v2) {
            return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
        },
        cross(v1, v2) {
            var a=v1[0]
            var b=v1[1]
            var c=v1[2]
            var d=v2[0]
            var e=v2[1]
            var f=v2[2]
            return [
                b * f - c * e,
                c * d - a * f,
                a * e - b * d,
            ]
        },
        x(v1, v) {
            return [
                v1[0] * v,
                v1[1] * v,
                v1[2] * v,
            ]
        },
        div(v1, v) {
            return [
                v1[0] / v,
                v1[1] / v,
                v1[2] / v,
            ]
        },
        sub(v1,v2) {
            return [
                v1[0] - v2[0],
                v1[1] - v2[1],
                v1[2] - v2[2],
            ]
        },
        add(v1, v2) {
            return [
                v1[0] + v2[0],
                v1[1] + v2[1],
                v1[2] + v2[2],
            ]
        },
        equal(v1,v2) {
            return v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2]
        },
        xMatrix(v1, v2) {
            return [
                v1[0]*v2[0][0]+v1[1]*v2[1][0]+v1[2]*v2[2][0]+v2[3][0],
                v1[0]*v2[0][1]+v1[1]*v2[1][1]+v1[2]*v2[2][1]+v2[3][1],
                v1[0]*v2[0][2]+v1[1]*v2[1][2]+v1[2]*v2[2][2]+v2[3][2],
                v1[0]*v2[0][3]+v1[1]*v2[1][3]+v1[2]*v2[2][3]+v2[3][3],
            ]
            var v1 = [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [v1[0], v1[1], v1[2], 1],
            ]
            var v = js3d.matrix16.x(v1,v2)
            return [v[3][0], v[3][1], v[3][2], v[3][3]]
        }
    },
    vector4: {
        clone(v) {
            return [
                v[0],
                v[1],
                v[2],
                v[3],
            ]
        }
    },
    vertex : {
        boxGeometryUnit() {
            var v8 = [
                [1, 1, 1],//0
                [1, 1, -1],//1
                [-1, 1, -1],//2
                [-1, 1, 1],//3
                [1, -1, 1],//4
                [1, -1, -1],//5
                [-1, -1, -1],//6
                [-1, -1, 1],//7
            ]
            var f4 = [
                [0, 1, 2, 3],//上
                [4, 7, 6, 5],//下
                [2, 1, 5, 6],//前
                [0, 3, 7, 4],//后
                [3, 2, 6, 7],//左
                [1, 0, 4, 5],//右
            ]
            var f3 = []
            for (var i=0; i<f4.length; ++i) {
                var f = f4[i]
                f3.push([
                    f[0],
                    f[1],
                    f[2],
                ])
                f3.push([
                    f[0],
                    f[2],
                    f[3],
                ])
            }
            return {
                v: v8,
                f4: f4,
                f3: f3,
            }
        },
        // radiusTop	可选。此属性定义圆柱体顶部圆半径。默认值是 20
        // radiusBottom	可选。此属性定义圆柱体底部圆半径。默认值是 20
        // height	可选。此属性定义圆柱体的高度。默认值是 100
        // radiusSegments	可选。此属性定义圆柱体的上下部的圆截面分成多少段。默认值是 8
        cylinderGeometry(radiusTop,height,radiusSegments) {
            if (radiusSegments<3) radiusSegments = 3
            var vertexArray = []
            vertexArray.push([0,height/2,0])
            var topBase = [radiusTop,height/2,0]
            var radiusStep = Math.PI/radiusSegments*2
            for (var i=0; i<radiusSegments; ++i) {
                var r = i*radiusStep
                var m = js3d.matrix16.buildRotateMatrixY(r)
                vertexArray.push(js3d.vector3.xMatrix(topBase,m))
            }
            var topVertexCnt = vertexArray.length
            for (var i=0; i<topVertexCnt; ++i) {
                vertexArray.push(js3d.vector3.add(vertexArray[i],[0,-height,0]))
            }
            var f3 = []
            var f4 = []
            var preI = topVertexCnt-1
            for (var i=1; i<topVertexCnt; ++i) {
                f3.push([i,0,preI])
                f3.push([preI+topVertexCnt,0,i+topVertexCnt])
                f4.push([i,preI,preI+topVertexCnt,i+topVertexCnt])
                preI = i
            }
            //var i = 9
            for (var i=0; i<f4.length; ++i) {
                var f = f4[i]
                f3.push([
                    f[0],
                    f[1],
                    f[2],
                ])
                f3.push([
                    f[0],
                    f[2],
                    f[3],
                 ])
            }
            return {
                v: vertexArray,
                f4: f4,
                f3: f3,
            }
        },
        boxGeometry(x, y, z, wx, wy, wz) {
            var m = this.boxGeometryUnit()
            for (var i = 0; i < m.v.length; ++i) {
                m.v[i][0] = m.v[i][0] * wx / 2 + x
                m.v[i][1] = m.v[i][1] * wy / 2 + y
                m.v[i][2] = m.v[i][2] * wz / 2 + z
            }
            return m
        },
    },
    color: {
        colorHexRGB2UnitRGBA(c) {
            return [
                ((c & 0xff0000) >> 16) / 0xff,
                ((c & 0x00ff00) >> 8) / 0xff,
                ((c & 0x0000ff)) / 0xff,
                1,
            ]
        },
    },
    matrix16 : {
        x(a, b) {
            var result = []
            for (var r = 0; r < 4; ++r) {
                var line = []
                for (var c = 0; c < 4; ++c) {
                var n = 0
                for (var c1 = 0; c1 < 4; ++c1) {
                    n += a[r][c1] * b[c1][c]
                }
                line.push(n)
                }
                result.push(line)
            }
            return result
        },
        buildUnitMatrix() {
            return [
                [1,0,0,0],
                [0,1,0,0],
                [0,0,1,0],
                [0,0,0,1],
            ]
        },
        buildAspectMatrix(aspect,fov,near,far) {
            var ver = 3
            if (ver==1) {
                var aspectI = 1/aspect
                var tanFovI = 1/Math.tan(fov/2)
                var disSubI = 1/(near-far)
                return [
                    [aspectI*tanFovI,0,0,0],
                    [0,tanFovI,0,0],
                    [0,0,(-far-near)*disSubI,2*near*far*disSubI],
                    [0,0,1,0],
                ]
            } else if (ver==2) {
                var aspectI = 1/aspect
                var cotFov = 1/Math.tan(fov/2)
                var disSubI = 1/(far-near)
                return [
                    [aspectI*cotFov,0,0,0],
                    [0,cotFov,0,0],
                    [0,0,far*disSubI,-near*far*disSubI],
                    [0,0,1,0],
                ]
            } else if (ver==3) {
                // NB: This creates 'uniform' perspective projection matrix,
                // which depth range [-1,1], right-handed rules
                //
                // [ A   0   C   0  ]
                // [ 0   B   D   0  ]
                // [ 0   0   q   qn ]
                // [ 0   0   -1  0  ]
                //
                // A  = 2 * near / (right - left)
                // B  = 2 * near / (top - bottom)
                // C  = (right + left) / (right - left)
                // D  = (top + bottom) / (top - bottom)
                // q  =     - (far + near) / (far - near)
                // qn = - 2 * (far * near) / (far - near)
                var halfFOV = fov*0.5;
                var tanHalfFOV = Math.tan(halfFOV);

                var top    = tanHalfFOV * near;
                var bottom = -top;
                var right  = top * aspect;
                var left   = -right;

                var width  = right - left;
                var height = top - bottom;

                var inv_w = 1 / width;
                var inv_h = 1 / height;
                var inv_d = 1 / (near - far);

                var A = 2 * near * inv_w;
                var B = 2 * near * inv_h;
                var C = 0//(right+left) * inv_w;
                var D = 0//(top+bottom) * inv_h;
                var q = (near==far) ? 1 : far*inv_d;
                var qn= near * q;

                return [
                    [A,0,C,0],
                    [0,B,D,0],
                    [0,0,q,qn],
                    [0,0,1,0],
                ]

                // {
                //     m_matProj.Zero();
                //     m_matProj._00 = A;
                //     m_matProj._11 = B;
                //     m_matProj._22 = q;
                //     m_matProj._32 = qn;
                //     m_matProj._23 = -1;
                //     //m_matProj._20 = C;
                //     //m_matProj._21 = D;
                // }
            }
        },
        buildScaleMatrix(x,y,z) {
            return [
                [x, 0, 0, 0],
                [0, y, 0, 0],
                [0, 0, z, 0],
                [0, 0, 0, 1],
            ]
        },
        buildTransformMatrix(x,y,z) {
            return [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [x, y, z, 1],
            ]
        },
        buildNormalize() {
            return [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
            ]
        },
        buildRotateMatrixX(radian) {
            var cos = Math.cos(radian)
            var sin = Math.sin(radian)
            return [
                [1, 0, 0, 0],
                [0, cos, -sin, 0],
                [0, sin, cos, 0],
                [0, 0, 0, 1],
            ]
        },
        buildRotateMatrixY(radian) {
            var cos = Math.cos(-radian)
            var sin = Math.sin(-radian)
            return [
                [cos, 0, sin, 0],
                [0, 1, 0, 0],
                [-sin, 0, cos, 0],
                [0, 0, 0, 1],
            ]
        },
        buildRotateMatrixZ(radian) {
            var cos = Math.cos(-radian)
            var sin = Math.sin(-radian)
            return [
                [cos, -sin, 0, 0],
                [sin, cos, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
            ]
        },
        buildLookAtMatrix(eye,lookat,up) {
            var ver = 2
            var eyex = eye[0]
            var eyey = eye[1]
            var eyez = eye[2]
            var centerx = lookat[0]
            var centery = lookat[1]
            var centerz = lookat[2]
            var upx = up[0]
            var upy = up[1]
            var upz = up[2]
            var m = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
            var x = [0, 0, 0]
            var y = [0, 0, 0]
            var z = [0, 0, 0]
        
            /* Difference eye and center vectors to make Z vector. */
            z[0] = centerx - eyex;
            z[1] = centery - eyey;
            z[2] = centerz - eyez;
            /* Normalize Z. */
            var mag = Math.sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
            if (mag) {
                z[0] /= mag;
                z[1] /= mag;
                z[2] /= mag;
            }
        
            /* Up vector makes Y vector. */
            y[0] = upx;
            y[1] = upy;
            y[2] = upz;
        
            /* X vector = Y cross Z. */
            x[0] = y[1] * z[2] - y[2] * z[1];
            x[1] = -y[0] * z[2] + y[2] * z[0];
            x[2] = y[0] * z[1] - y[1] * z[0];
        
            /* Recompute Y = Z cross X. */
            y[0] = z[1] * x[2] - z[2] * x[1];
            y[1] = -z[0] * x[2] + z[2] * x[0];
            y[2] = z[0] * x[1] - z[1] * x[0];
        
            /* Normalize X. */
            mag = Math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
            if (mag) {
                x[0] /= mag;
                x[1] /= mag;
                x[2] /= mag;
            }
        
            /* Normalize Y. */
            mag = Math.sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
            if (mag) {
                y[0] /= mag;
                y[1] /= mag;
                y[2] /= mag;
            }
        
            if (ver==1) {
                /* Build resulting view matrix. 行主序*/
                m[0 * 4 + 0] = x[0]; m[0 * 4 + 1] = x[1];
                m[0 * 4 + 2] = x[2]; m[0 * 4 + 3] = -x[0] * eyex + -x[1] * eyey + -x[2] * eyez;
            
                m[1 * 4 + 0] = y[0]; m[1 * 4 + 1] = y[1];
                m[1 * 4 + 2] = y[2]; m[1 * 4 + 3] = -y[0] * eyex + -y[1] * eyey + -y[2] * eyez;
            
                m[2 * 4 + 0] = z[0]; m[2 * 4 + 1] = z[1];
                m[2 * 4 + 2] = z[2]; m[2 * 4 + 3] = -z[0] * eyex + -z[1] * eyey + -z[2] * eyez;
            
                m[3 * 4 + 0] = 0.0; m[3 * 4 + 1] = 0.0; m[3 * 4 + 2] = 0.0; m[3 * 4 + 3] = 1.0;
    
                if (false) {
                    return [
                        [m[0], m[4], m[8], m[12]],
                        [m[1], m[5], m[9], m[13]],
                        [m[2], m[6], m[10], m[14]],
                        [m[3], m[7], m[11], m[15]],
                    ]
                } else {
                    return [
                        [m[0], m[1], m[2], m[3]],
                        [m[4], m[5], m[6], m[7]],
                        [m[8], m[9], m[10], m[11]],
                        [m[12], m[13], m[14], m[15]],
                    ]
                }
            } else if (ver==2) {
                var T = this.buildTransformMatrix(-eyex,-eyey,-eyez)//平移到摄像机坐标系下
                var R = [
                    [x[0],x[1],x[2],0],
                    [y[0],y[1],y[2],0],
                    [z[0],z[1],z[2],0],
                    [0,0,0,1],
                ]
                return this.x(T,R)
            }
        }
    },
    intersect: {
        // 空间点 sp 起点 sq终点
        // 三角形空间点 a b c 
        // 输出参数 t 
        line2triangle(sp, sq, a, b, c)
        {
            var ab = js3d.vector3.sub(b,a)
            var ac = js3d.vector3.sub(c,a)
            // 终点->起点 的向量, 不是 起点->终点
            var qp = js3d.vector3.sub(sp, sq);
    
            // Compute triangle normal. Can be precalculated or cached if
            // intersecting multiple segments against the same triangle
            // 法线方向为三角形正面方向 
            var norm = js3d.vector3.cross(ab, ac);
    
            // Compute denominator d. If d <= 0, segment is parallel to or points
            // away from triangle, so exit early
            // d = 0.0 说明 qp 和 norm 垂直，说明三角形和 qp 平行。
            // d < 0.0 说明 qp 和 norm 是钝角 说明是从三角形的背面 进入和三角形相交的 
            var d = js3d.vector3.dot(qp, norm);
            if (d <= 0) return null;
    
            // Compute intersection t value of pq with plane of triangle. A ray
            // intersects if 0 <= t. Segment intersects if 0 <= t <= 1. Delay
            // dividing by d until intersection has been found to pierce triangle
            var ap = js3d.vector3.sub(sp, a);
            var t = js3d.vector3.dot(ap, norm);
            if (t < 0) return null;
            if (t > d) return null; // For segment; exclude this code line for a ray test 
    
            // Compute barycentric coordinate components and test if within bounds
            var e = js3d.vector3.cross(qp, ap);
            var v = js3d.vector3.dot(ac, e);
            if (v < 0 || v > d) return null;
            var w = -js3d.vector3.dot(ab, e);
            if (w < 0 || v + w > d) return null;
    
            // Segment/ray intersects triangle. Perform delayed division
            t /= d;
    
            return t;
        }
    },
    device: {
        context:null,
        imageData:null,
        imageSize:[0,0],
        Flush:function(data,width,height,targetX,targetY,targetWidth,targetHeight) {
            if (this.context==null) return
            if (this.imageData==null || this.imageSize[0]!=width || this.imageSize[1]!=height) {                
                this.imageData = this.context.createImageData(width,height)
            }
            var imageData = this.imageData
            var src = data
            var dist = imageData.data
            var i = 0;
            var j = 0;
            for (; i<src.length; i+=5,j+=4) {
                dist[j+0] = 255*src[i+0]
                dist[j+1] = 255*src[i+1]
                dist[j+2] = 255*src[i+2]
                dist[j+3] = 255*src[i+3]
            }
            this.FlashImageData(imageData,targetX,targetY,targetWidth,targetHeight)
        },
        FlashImageData(imageData,targetX,targetY,targetWidth,targetHeight) {
            if (this.context==null) return
            this.context.putImageData(imageData,0,0,targetX,targetY,targetWidth,targetHeight)
        }
    },
    engine: {
        device: null,
        size: [0,0],
        far:10000,
        near:2,
        espectMatrix : null,
        dataLine: null,
        matrixStack:[],
        eye : [0,100,-200],
        lookat : [0,100,0],
        up : [0,1,0],
        viewMatrix: null,
        viewMatrixEspect: null,
        Init: function(width,height,device) {
            this.device = device
            this.size[0] = width
            this.size[1] = height
            this.dataLine = []
            for (var i=0; i<width*height; ++i) {
                this.dataLine.push(0)
                this.dataLine.push(0)
                this.dataLine.push(0)
                this.dataLine.push(0)
                this.dataLine.push(null)
            }
            this.espectMatrix = js3d.matrix16.buildAspectMatrix(this.size[0]/this.size[1],(Math.PI/3),this.near,this.far)
            this.CleanBuffer()
        },
        ClipDEx(d1,d2,a1,a2,a3,near,far) {
            if (d1.v[a1]<near) {
                if (d2.v[a1]<0) return false//丢弃线段
                d1 = {
                    v:js3d.vector3.clone(d1.v),
                    c:js3d.vector4.clone(d1.c),
                }
                d2 = {
                    v:js3d.vector3.clone(d2.v),
                    c:js3d.vector4.clone(d2.c),
                }
                if (d2.v[a1] - d1.v[a1]==0) return null//丢弃三角形
                var t = (near-d1.v[a1]) / (d2.v[a1] - d1.v[a1])
                d1.v[a1] = near
                d1.v[a2] = d1.v[a2] + t * (d2.v[a2]-d1.v[a2])
                d1.v[a3] = d1.v[a3] + t * (d2.v[a3]-d1.v[a3])
                return d1//d1已经修改
            } else if (d1.v[1]>=far) {
                if (d2.v[a1]>=far) return false//丢弃线段
                d1 = {
                    v:js3d.vector3.clone(d1.v),
                    c:js3d.vector4.clone(d1.c),
                }
                d2 = {
                    v:js3d.vector3.clone(d2.v),
                    c:js3d.vector4.clone(d2.c),
                }
                if (d1.v[a1] - d2.v[a1]==0) return null//丢弃三角形
                var t = (d1.v[a1]-far) / (d1.v[a1] - d2.v[a1])
                d1.v[a1] = far-1
                d1.v[a2] = d1.v[a2] + t * (d2.v[a2]-d1.v[a2])
                d1.v[a3] = d1.v[a3] + t * (d2.v[a3]-d1.v[a3])
                return d1
            } else {
                return true//没有变化
            }
        },
        ClipD(d1,d2) {
            if (d1.v[1]<0) {
                if (d2.v[1]<0) return false
                d1 = {
                    v:js3d.vector3.clone(d1.v),
                    c:js3d.vector4.clone(d1.c),
                }
                d2 = {
                    v:js3d.vector3.clone(d2.v),
                    c:js3d.vector4.clone(d2.c),
                }
                if (d2.v[1] - d1.v[1]==0) return false
                var t = -d1.v[1] / (d2.v[1] - d1.v[1])
                d1.v[1] = 0
                d1.v[0] = d1.v[0] + t * (d2.v[0]-d1.v[0])
                d1.v[2] = d1.v[2] + t * (d2.v[2]-d1.v[2])
            } else if (d1.v[1]>=this.size[1]) {
                if (d2.v[1]>=this.size[1]) return false
                d1 = {
                    v:js3d.vector3.clone(d1.v),
                    c:js3d.vector4.clone(d1.c),
                }
                d2 = {
                    v:js3d.vector3.clone(d2.v),
                    c:js3d.vector4.clone(d2.c),
                }
                var t = (d1.v[1]-this.size[1]) / (d1.v[1] - d2.v[1])
                d1.v[1] = this.size[1]-1
                d1.v[0] = d1.v[0] + t * (d2.v[0]-d1.v[0])
                d1.v[2] = d1.v[2] + t * (d2.v[2]-d1.v[2])
            }
            return [d1,d2]
        },
        CalcScreenLine(d1,d2,data) {
            var _d0 = Math.abs(d1.v[0]-d2.v[0])
            var _d1 = Math.abs(d1.v[1]-d2.v[1])
            if (_d0==0 && _d1==0) return
            var d = this.ClipD(d1,d2)
            if (!d) return
            d1 = d[0]; d2 = d[1];
            var d = this.ClipD(d2,d1)
            if (!d) return
            d1 = d[1]; d2 = d[0];
            var z1 = d1
            var z2 = d2
            var v1 = d1.v
            var v2 = d2.v
            var c1 = d1.c
            var c2 = d2.c
            var d0 = v2[0]-v1[0]
            var d1 = v2[1]-v1[1]
            var d2 = v2[2]-v1[2]
            var _d0 = Math.abs(d0)
            var _d1 = Math.abs(d1)
            if (_d0==0 && _d1==0) return
            ++ref.line
            var t0 = 0
            var t1 = 0
            var t2 = 0
            var tc = [0,0,0,0]
            var times = 0
            if (_d0>_d1) {
                t0 = d0 / _d0
                t1 = d1 / _d0
                t2 = d2 / _d0
                for (var i=0; i<4; ++i)
                    tc[i] = (c2[i]-c1[i]) / _d0
                times = _d0
            } else {
                t0 = d0 / _d1
                t1 = d1 / _d1
                t2 = d2 / _d1
                for (var i=0; i<4; ++i)
                    tc[i] = (c2[i]-c1[i]) / _d1
                times = _d1
            }
            var x = v1[0]
            var y = v1[1]
            var z = v1[2]
            var c = [
                c1[0],
                c1[1],
                c1[2],
                c1[3],
            ]
            while(times>0){
                var w = Math.floor(x)
                var h = Math.floor(y)
                var hs = h+''
                var d = data[hs]
                if (d==null) {
                    d = {
                        h:h,
                        w:[{
                            w:w,
                            z:z,
                            c:[
                                c[0],
                                c[1],
                                c[2],
                                c[3],
                            ],
                        }],
                    }
                    data[hs] = d
                }
                ++ref.linePixel
                d.w.push({
                    w:w,
                    z:z,
                    c:[
                        c[0],
                        c[1],
                        c[2],
                        c[3],
                    ],
                })
                x += t0
                y += t1
                z += t2
                for (var j=0; j<4; ++j)
                    c[j] += tc[j]
                --times
            }
        },
        DrawScreenTriangles(vertexArray,lightStrong) {
            var range = [this.size[0],this.size[1],1]
            for (var i=0; i<3; ++i) {
                if (vertexArray[0].v[i]<0 && vertexArray[1].v[i]<0 && vertexArray[2].v[i]<0) return
                if (vertexArray[0].v[i]>=range[i] && vertexArray[1].v[i]>=range[i] && vertexArray[2].v[i]>=range[i]) return
            }
            var last = null
            var data = {}
            ++ref.trangle
            for (var i=0; i<vertexArray.length; ++i) {
                if (last) {
                    this.CalcScreenLine(last,vertexArray[i],data)
                }
                last = vertexArray[i]
            }
            this.CalcScreenLine(last,vertexArray[0],data)
            for (var k in data) {
                var d = data[k]
                var h = d.h
                var warray = d.w
                warray.sort(function(v1,v2){
                    if (v1.w<v2.w) return -1
                    if (v1.w==v2.w) return 0
                    return 1
                })
                var lastw = null
                for (var i in warray) {
                    var cw = d.w[i]
                    if (lastw!=null && (lastw.w+1)<=cw.w) {
                        var tc = [
                            (cw.c[0]-lastw.c[0])/(cw.w-lastw.w),
                            (cw.c[1]-lastw.c[1])/(cw.w-lastw.w),
                            (cw.c[2]-lastw.c[2])/(cw.w-lastw.w),
                            (cw.c[3]-lastw.c[3])/(cw.w-lastw.w),
                        ]
                        var c = [
                            lastw.c[0],
                            lastw.c[1],
                            lastw.c[2],
                            lastw.c[3],
                        ]
                        var dw = 1
                        var dz = (cw.z-lastw.z)/(cw.w-lastw.w)
                        var w = lastw.w
                        var z = lastw.z
                        while(true) {
                            w += dw
                            z += dz
                            if (w>cw.w)
                                break
                            for (var j=0; j<4; ++j)
                                c[j] += tc[j]
                            this.SetPixel(w,h,z,c[0],c[1],c[2],c[3],lightStrong)
                        }
                    } else {
                        var w = cw.w
                        var c = cw.c
                        var z = cw.z
                        this.SetPixel(w,h,z,c[0],c[1],c[2],c[3],lightStrong)
                    }
                    lastw = cw
                }
            }
        },
        SetPixel(w,h,z,r,g,b,a,lightStrong) {
            if (lightStrong<0.3) lightStrong = 0.3
            h = this.size[1]-h
            if (h<0 || h>= this.size[1]) return
            if (w<0 || w>= this.size[0]) return
            var i=(h*this.size[0]+w)*5
            var src = this.dataLine
            if (src[i+4]!=null && src[i+4]<=z) return
            src[i+0] = r*lightStrong
            src[i+1] = g*lightStrong
            src[i+2] = b*lightStrong
            src[i+3] = a
            src[i+4] = z
        },
        DrawTriangles(vertexArray) {
            var m = this.GetLastMatrix()
            var finalVertexArray = []
            for (var i=0; i<vertexArray.length; ++i) {
                var d = vertexArray[i]
                var c = [d[3],d[4],d[5],d[6]]
                finalVertexArray.push({
                    v:d,
                    c:c,
                })
            }
            var cubeViewVertexArray = []
            var viewMatrix = this.viewMatrix
            if (m)
                viewMatrix = js3d.matrix16.x(m,viewMatrix)
            for (var i=0; i<finalVertexArray.length; ++i) {
                var d = finalVertexArray[i]
                var c = d.c
                var v = d.v
                v = js3d.vector3.xMatrix(v,viewMatrix)
                cubeViewVertexArray.push({
                    v:v,
                    c:c,
                })
            }
            //剪切到[near,far]
            var newVertexArray = []
            for (var i=0; i<cubeViewVertexArray.length; ++i) {
                var j=i-1
                var k=i+1
                if (j<0) j=cubeViewVertexArray.length-1
                if (k>=cubeViewVertexArray.length) k=0
                var r = this.ClipDEx(cubeViewVertexArray[i],cubeViewVertexArray[j],2,0,1,this.near,this.far*100)
                if (r===true) {
                    //没有变化
                    newVertexArray.push(cubeViewVertexArray[i])
                } else {
                    if (r===false) {
                        //丢弃线段
                    } else {
                        newVertexArray.push(r)
                    }
                    var r = this.ClipDEx(cubeViewVertexArray[i],cubeViewVertexArray[k],2,0,1,this.near,this.far*100)
                    if (r===false) {
                        //丢弃线段
                    } else if (r===true) {
                        //没有变化//错误情况
                        console.log('error-----------')
                    } else {
                        newVertexArray.push(r)
                    }
                }
            }
            cubeViewVertexArray = newVertexArray
            if (cubeViewVertexArray.length<3) return
            //投影
            var screenVertexArray = []
            for (var i=0; i<cubeViewVertexArray.length; ++i) {
                var d = cubeViewVertexArray[i]
                var c = d.c
                var v = d.v
                v = js3d.vector3.xMatrix(v,this.espectMatrix)
                if (v[3]>0&&v[3]<0.0001) v[3]=0.0001
                if (v[3]<0&&v[3]>-0.0001) v[3]=-0.0001
                var wAbs = Math.abs(v[3])
                if (wAbs>0) {
                    v[0] /= wAbs
                    v[1] /= wAbs
                    v[2] /= -wAbs
                }
                v[0] = this.size[0]*(v[0]+0.5)
                v[1] = this.size[1]*(v[1]+0.5)
                screenVertexArray.push({
                    v:v,
                    c:c,
                })
            }
            //剔除背面
            var v0 = screenVertexArray[0].v
            var v1 = screenVertexArray[1].v
            var v2 = screenVertexArray[2].v
            var d1 = js3d.vector2.normalize(js3d.vector2.sub(v1,v0))
            var d2 = js3d.vector2.normalize(js3d.vector2.sub(v2,v0))
            var faceN = js3d.vector2.cross(d1,d2)
            if (faceN>=0) return
            //需要渲染//开始计算光照
            if (m) {
                for (var i=0; i<finalVertexArray.length; ++i) {
                    finalVertexArray[i].v = js3d.vector3.xMatrix(finalVertexArray[i].v,m)
                }
            }
            //计算法线
            var v0 = finalVertexArray[0].v
            var v1 = finalVertexArray[1].v
            var v2 = finalVertexArray[2].v
            var d1 = js3d.vector3.normalize(js3d.vector3.sub(v1,v0))
            var d2 = js3d.vector3.normalize(js3d.vector3.sub(v2,v0))
            var faceN = js3d.vector3.cross(d1,d2)
            faceN = js3d.vector3.normalize(faceN)
            //计算平行光强度
            var light = js3d.vector3.normalize([-3,-2,3])
            var dot = js3d.vector3.dot(light,faceN)
            var lightStrong = Math.max(0,-dot*2)
            for (var i=2; i<screenVertexArray.length; ++i) {
                var v = [
                    screenVertexArray[0],
                    screenVertexArray[i-1],
                    screenVertexArray[i]
                ]
                this.DrawScreenTriangles(v,lightStrong)
            }
        },
        CleanBuffer() {
            var c = js3d.color.colorHexRGB2UnitRGBA(0xf7d9aa)
            var d = this.dataLine
            var len = d.length
            for (var i=0; i<len; i+=5) {
                d[i+0] = c[0]
                d[i+1] = c[1]
                d[i+2] = c[2]
                d[i+3] = 1
                d[i+4] = null
            }
        },
        Clean: function() {
            this.CleanBuffer()
            this.ClearMatrix()
        },
        Flush: function() {
            if (this.device==null) return
            this.device.Flush(this.dataLine,this.size[0],this.size[1],0,0,this.size[0],this.size[1])
        },
        PushMatrix: function(m) {
            if (!m) return
            var lastIdx = this.matrixStack.length-1;
            var cur = this.matrixStack[lastIdx]
            if (cur == null)
                cur = js3d.matrix16.buildNormalize()
            this.matrixStack.push(js3d.matrix16.x(cur,m))
        },
        GetLastMatrix: function() {
            var lastIdx = this.matrixStack.length-1;
            return this.matrixStack[lastIdx]
        },
        ClearMatrix() {
            this.matrixStack = []
        },
        SetMatrix(m) {
            this.ClearMatrix()
            this.PushMatrix(m)
        },
        BeforeRender() {
            //视锥
            var eye = this.eye
            var lookat = this.lookat
            var up = this.up
            this.viewMatrix = js3d.matrix16.buildLookAtMatrix(eye,lookat,up)
            this.viewMatrixEspect = js3d.matrix16.x(this.viewMatrix,this.espectMatrix)
        }
    }
}
