module bch(
	input clk,
	input rstn,
	input mode,
	input [1:0] code,module bch(
	input clk,
	input rstn,
	input mode,
	input [1:0] code,
	input set,
	input [63:0] idata,
	output ready,
	output finish,
	output [9:0] odata
);
// -------------------------------------------------------------------------------------------------------------------------- //
// Declaration
// -------------------------------------------------------------------------------------------------------------------------- //
	reg       ready_reg;
	assign    ready = ready_reg;

	reg [1:0] code_reg; // 1: (63,51); 2: (255,239); 3: (1023,983)
	reg       mode_reg; // 0: hard-decisioan decoding; 1: soft-decision decoding

// counter
	reg [2:0] cnt, cnt_n;

// LLR
	reg [1:1023] llr, llr_n; // 1: MSB ;  1023: LSB
	reg [7:0]    idata_in;

// CRC
	reg [0:4] w_CRC, r_CRC; // 0: MSB ; 4: LSB
	// reg [0:4] 
	reg [2:0] w_CRC_remainder_0, w_CRC_remainder_1, w_CRC_remainder_2, w_CRC_remainder_3, w_CRC_remainder_4, w_CRC_remainder_5, w_CRC_remainder_6, w_CRC_remainder_7, w_CRC_remainder_8, w_CRC_remainder_9, w_CRC_remainder_10;
    
// BM
	reg [5:0] syndromes [0:3], syndromes_n [0:3];
	reg [5:0] sigma [0:2], sigma_n [0:2];
	reg [5:0] sigma_tmp [0:2], sigma_n_tmp [0:2];
	reg [5:0] aux_sigma [0:2], aux_sigma_n [0:2];

	reg [5:0] delta, delta_n;
	reg [2:0] idx_cnt, idx_cnt_n;
	reg [5:0] gamma, gamma_n;
	reg [2:0] k, k_n;

	wire [2:0] t;
	assign t = mode == 2'd3 ? 4 : 2;
	wire [2:0] r1, r2, base_idx;
	assign r1 = (t<<1) - idx_cnt;
	assign r2 = r1+t > (t<<1) ? (t<<1) : r1+t;
	assign base_idx = (t<<1)-r2;
	wire [2:0] num;
	assign num = r2-r1;
	
	reg [5:0] delta_poly, delta_poly_n; 
	wire [5:0] syndrome_poly, syndrome_power, sigma_poly, sigma_power, delta_poly_w;
	wire [6:0] overflow_power;
	assign overflow_power = (syndrome_power + sigma_power) >= 63 ? syndrome_power + sigma_power - 63 : syndrome_power + sigma_power;
	assign syndrome_poly = syndromes[idx_cnt-cnt];
	assign sigma_poly = sigma[cnt];
	assign finish = 0;
	// assign ready = 1;
	// assign finish = 0;
	reg [2:0] c_state, n_state;
	localparam S_IDLE  = 3'd0,
			   S_FETCH = 3'd1,
			   S_CRC   = 3'd2,
			   S_BM_delta = 3'd3,
			   S_BM_update_1 = 3'd4,
			   S_BM_update_2 = 3'd5,
			   S_CHIEN =  3'd6,
			   S_END   = 3'd7;


// BM algorithm
// -------------------------------------------------------------------------------------------------------------------------- //
	LUT lut_unit(.code(code), .poly(c_state == S_BM_delta ? syndrome_poly : c_state == S_BM_update_1 ? gamma : delta_poly), .power(syndrome_power));
	LUT sigma_lut_unit(.code(code), .poly(c_state == S_BM_update_2 ? aux_sigma[cnt] : sigma_poly), .power(sigma_power));
	LUT_inv delta_lut_unit(.code(code), .power(c_state == S_BM_delta ? delta_n : overflow_power), .poly(delta_poly_w));

	integer i = 0;
	always @(*) begin
		idx_cnt_n = idx_cnt;
		delta_n = delta;
		delta_poly_n = delta_poly;
		gamma_n = gamma;
		k_n = k;
		for(i=0; i<4; i++) syndromes_n[i] = syndromes[i];
		for(i=0; i<3; i++) begin
			sigma_n[i] = sigma[i];
			aux_sigma_n[i] = aux_sigma[i];
		end

		case(c_state)
			S_BM_delta: begin // 算 delta
				idx_cnt_n = cnt == num ? idx_cnt == (t<<1) ? idx_cnt : idx_cnt + 1 : idx_cnt; // 記現在是在處理第幾個symptome
				delta_n = sigma_poly==0 || syndrome_poly==0 ? -1 : overflow_power; // 查表加法 （要考慮是否為0，不然會算錯), overflow要特別處理，因為overflow會自動 mod 64 但是我們要的是 mod 63
				delta_poly_n = cnt==0 ? delta_poly_w : delta_poly_w ^ delta_poly; // 每個 cycle 都跟前一個 cycle 做完的 delta 做 XOR ，如果 cnt==0 就不做
				for(i=0; i<3; i++) sigma_n_tmp[i] = sigma[i]; // 先記得 sigma 函數，後面有可能要指派給 aux_sigma
			end
			S_BM_update_1: begin  // 先做 sigma_n = gamma * sigma + delta * aux_sigma 中的前一項 gamma * sigma (t+1個cycle，因為sigma有t+1個根)
				sigma_n[cnt] = sigma[cnt] == 0 || gamma == 0 ? 0 : delta_poly_w;
			end
			S_BM_update_2: begin // 再做 sigma_n = gamma * sigma + delta * aux_sigma 中的後一項 delta * aux_sigma (t+1個cycle，因為sigma有t+1個根)
				sigma_n[cnt] = aux_sigma[cnt] == 0 || delta_poly == 0 ? sigma[cnt] : sigma[cnt] + delta_poly_w;
				if(delta_poly != 0 && cnt==t+1 && k >= 0) begin
					for(i=0; i<3; i++) begin
						aux_sigma_n[i] = i == 0 ? 0 : sigma_tmp[i-1];
					end
					gamma_n = delta_poly;
					k_n = -k-1;
				end
				else if (cnt==t+1) begin
					for(i=0; i<3; i++) begin
						aux_sigma_n[i] = i == 0 ? 0 : aux_sigma[i-1];
					end
					k_n = k+1;
				end
			end
			S_CHIEN: begin
				$display("=================FINISH BM==================");
				$write("sigma: ");
				for(i=0; i<3; i++) begin
					$write(sigma[i], " ");
				end
				$write("\n");
				$display("=================FINISH BM==================");
			end

		endcase
	end

	always @(posedge clk) begin
		if(~rstn) begin
			idx_cnt <= 0;
			delta <= 0;
			delta_poly <= 0;
			gamma <= 1;
			k <= 0;
			for(i=0; i<3; i++) begin
				sigma[i] <= 0;
				aux_sigma[i] <= 0;
				sigma_tmp[i] <= 0;
			end
			sigma[0] <= 1;
			sigma_tmp[0] <= 0;
			aux_sigma[1] <= 1;

			// From tb p100-1
			syndromes[0] <= 6'd2;
			syndromes[1] <= 6'd4;
			syndromes[2] <= 6'd32;
			syndromes[3] <= 6'd16;
			
			// From tb p100-2
			// syndromes[0] <= 6'd31;
			// syndromes[1] <= 6'd26;
			// syndromes[2] <= 6'd13;
			// syndromes[3] <= 6'd11;
		end
		else begin
			idx_cnt <= idx_cnt_n;
			delta <= delta_n;
			delta_poly <= delta_poly_n;
			gamma <= gamma_n;
			k <= k_n;
			for(i=0; i<3; i++) begin
				sigma[i] <= sigma_n[i];
				aux_sigma[i] <= aux_sigma_n[i];
				sigma_tmp[i] <= sigma_n_tmp[i];
			end

			for(i=0; i<4; i++) begin
				syndromes[i] <= syndromes_n[i];
			end
		end
	end
// -------------------------------------------------------------------------------------------------------------------------- //
// Logic
// -------------------------------------------------------------------------------------------------------------------------- //
// FSM
	always @(posedge clk) begin
		if(~rstn) begin
			c_state <= S_IDLE;
		end
		else begin
			c_state <= n_state;
		end
	end
	always @* begin
		n_state = c_state;
		case (c_state)
			S_IDLE: n_state = S_FETCH;
			S_FETCH: n_state = (cnt == 7)? S_CRC:S_FETCH;
			S_CRC  : n_state = (cnt == 7)? S_BM_delta:S_CRC;
			S_BM_delta : n_state = cnt==num ? S_BM_update_1: S_BM_delta;
			S_BM_update_1 : n_state = cnt == t+1 ? S_BM_update_2 : S_BM_update_1;
			S_BM_update_2 : n_state = cnt == t+1 ? idx_cnt == (t<<1) ? S_CHIEN : S_BM_delta : S_BM_update_2;
			S_CHIEN: n_state = S_END;
			//default: 
		endcase
	end

// code info.
	always @(posedge clk) begin
		if(~rstn) begin
			code_reg <= 0;
			mode_reg <= 0;
		end
		else if(set) begin
			code_reg <= code;
			mode_reg <= mode;
		end
	end

// counter
	always @(posedge clk) begin
		if(~rstn) begin
			cnt <= 0;
		end
		else begin
			cnt <= cnt_n;
		end
	end

	always @* begin
		cnt_n = cnt;
		case(c_state)
			S_FETCH: cnt_n = (cnt == 7)? 0:cnt + 1;
			S_CRC: cnt_n = (cnt == 7)? 0:cnt + 1;
			S_BM_delta:  cnt_n = cnt == num ? 0 : cnt + 1;
			S_BM_update_1: cnt_n = (cnt == t+1) ? 0 : cnt + 1;
			S_BM_update_2: cnt_n = (cnt == t+1) ? 0 : cnt + 1;
			S_CHIEN: cnt_n = 0;
			S_END: cnt_n = 0;
		endcase
	end

// llr
	always @(posedge clk) begin
		if(~rstn) begin
			llr <= 0;
		end
		else begin
			llr <= llr_n;
		end
	end
	always @* begin
		llr_n = llr;
		case(c_state)
			S_FETCH: begin
				case(cnt)
					0: llr_n[1:7] = {idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					1: llr_n[8:15] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					2: llr_n[16:23] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					3: llr_n[24:31] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					4: llr_n[32:39] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					5: llr_n[40:47] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					6: llr_n[48:55] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					7: llr_n[56:63] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
				endcase
			end
		endcase
	end




// output 
	always @* begin
		if(c_state == S_FETCH) ready_reg = 1;
		else ready_reg = 0;
	end

endmodule

module LUT(
	input [1:0] code,
	input [5:0] poly,
	output reg [5:0] power
);
	always @* begin
		power = 6'd0;
		case(poly)
			6'd0: power = -6'd1;
			6'd1: power = 6'd0;
			6'd2: power = 6'd1;
			6'd4: power = 6'd2;
			6'd8: power = 6'd3;
			6'd16: power = 6'd4;
			6'd32: power = 6'd5;
			6'd3: power = 6'd6;
			6'd6: power = 6'd7;
			6'd12: power = 6'd8;
			6'd24: power = 6'd9;
			6'd48: power = 6'd10;
			6'd35: power = 6'd11;
			6'd5: power = 6'd12;
			6'd10: power = 6'd13;
			6'd20: power = 6'd14;
			6'd40: power = 6'd15;
			6'd19: power = 6'd16;
			6'd38: power = 6'd17;
			6'd15: power = 6'd18;
			6'd30: power = 6'd19;
			6'd60: power = 6'd20;
			6'd59: power = 6'd21;
			6'd53: power = 6'd22;
			6'd41: power = 6'd23;
			6'd17: power = 6'd24;
			6'd34: power = 6'd25;
			6'd7: power = 6'd26;
			6'd14: power = 6'd27;
			6'd28: power = 6'd28;
			6'd56: power = 6'd29;
			6'd51: power = 6'd30;
			6'd37: power = 6'd31;
			6'd9: power = 6'd32;
			6'd18: power = 6'd33;
			6'd36: power = 6'd34;
			6'd11: power = 6'd35;
			6'd22: power = 6'd36;
			6'd44: power = 6'd37;
			6'd27: power = 6'd38;
			6'd54: power = 6'd39;
			6'd47: power = 6'd40;
			6'd29: power = 6'd41;
			6'd58: power = 6'd42;
			6'd55: power = 6'd43;
			6'd45: power = 6'd44;
			6'd25: power = 6'd45;
			6'd50: power = 6'd46;
			6'd39: power = 6'd47;
			6'd13: power = 6'd48;
			6'd26: power = 6'd49;
			6'd52: power = 6'd50;
			6'd43: power = 6'd51;
			6'd21: power = 6'd52;
			6'd42: power = 6'd53;
			6'd23: power = 6'd54;
			6'd46: power = 6'd55;
			6'd31: power = 6'd56;
			6'd62: power = 6'd57;
			6'd63: power = 6'd58;
			6'd61: power = 6'd59;
			6'd57: power = 6'd60;
			6'd49: power = 6'd61;
			6'd33: power = 6'd62;
		endcase
	end
endmodule

module LUT_inv (
	input [1:0] code,
	input [5:0] power,
	output reg [5:0] poly
);
	always @(*) begin
		case(power)
			-6'd1: poly = 6'd0;
			6'd0: poly = 6'd1;
			6'd1: poly = 6'd2;
			6'd2: poly = 6'd4;
			6'd3: poly = 6'd8;
			6'd4: poly = 6'd16;
			6'd5: poly = 6'd32;
			6'd6: poly = 6'd3;
			6'd7: poly = 6'd6;
			6'd8: poly = 6'd12;
			6'd9: poly = 6'd24;
			6'd10: poly = 6'd48;
			6'd11: poly = 6'd35;
			6'd12: poly = 6'd5;
			6'd13: poly = 6'd10;
			6'd14: poly = 6'd20;
			6'd15: poly = 6'd40;
			6'd16: poly = 6'd19;
			6'd17: poly = 6'd38;
			6'd18: poly = 6'd15;
			6'd19: poly = 6'd30;
			6'd20: poly = 6'd60;
			6'd21: poly = 6'd59;
			6'd22: poly = 6'd53;
			6'd23: poly = 6'd41;
			6'd24: poly = 6'd17;
			6'd25: poly = 6'd34;
			6'd26: poly = 6'd7;
			6'd27: poly = 6'd14;
			6'd28: poly = 6'd28;
			6'd29: poly = 6'd56;
			6'd30: poly = 6'd51;
			6'd31: poly = 6'd37;
			6'd32: poly = 6'd9;
			6'd33: poly = 6'd18;
			6'd34: poly = 6'd36;
			6'd35: poly = 6'd11;
			6'd36: poly = 6'd22;
			6'd37: poly = 6'd44;
			6'd38: poly = 6'd27;
			6'd39: poly = 6'd54;
			6'd40: poly = 6'd47;
			6'd41: poly = 6'd29;
			6'd42: poly = 6'd58;
			6'd43: poly = 6'd55;
			6'd44: poly = 6'd45;
			6'd45: poly = 6'd25;
			6'd46: poly = 6'd50;
			6'd47: poly = 6'd39;
			6'd48: poly = 6'd13;
			6'd49: poly = 6'd26;
			6'd50: poly = 6'd52;
			6'd51: poly = 6'd43;
			6'd52: poly = 6'd21;
			6'd53: poly = 6'd42;
			6'd54: poly = 6'd23;
			6'd55: poly = 6'd46;
			6'd56: poly = 6'd31;
			6'd57: poly = 6'd62;
			6'd58: poly = 6'd63;
			6'd59: poly = 6'd61;
			6'd60: poly = 6'd57;
			6'd61: poly = 6'd49;
			6'd62: poly = 6'd33;
			default: poly = 6'd0;
		endcase
	end
endmodule
	input set,
	input [63:0] idata,
	output ready,
	output finish,
	output [9:0] odata
);
// -------------------------------------------------------------------------------------------------------------------------- //
// Declaration
// -------------------------------------------------------------------------------------------------------------------------- //
	reg       ready_reg;
	assign    ready = ready_reg;

	reg [1:0] code_reg; // 1: (63,51); 2: (255,239); 3: (1023,983)
	reg       mode_reg; // 0: hard-decision decoding; 1: soft-decision decoding

// counter
	reg [2:0] cnt, cnt_n;

// LLR
	reg [1:1023] llr, llr_n; // 1: MSB ;  1023: LSB
	reg [7:0]    idata_in;

// CRC
	reg [0:4] w_CRC, r_CRC; // 0: MSB ; 4: LSB
	// reg [0:4] 
	reg [2:0] w_CRC_remainder_0, w_CRC_remainder_1, w_CRC_remainder_2, w_CRC_remainder_3, w_CRC_remainder_4, w_CRC_remainder_5, w_CRC_remainder_6, w_CRC_remainder_7, w_CRC_remainder_8, w_CRC_remainder_9, w_CRC_remainder_10;
    

// assign finish = (cnt == 7 || cnt == 8)? 1:0;

	reg [2:0] c_state, n_state;
	localparam S_IDLE  = 3'b000,
			   S_FETCH = 3'b001,
			   S_CRC   = 3'b011,
			   S_END   = 3'b111;


// -------------------------------------------------------------------------------------------------------------------------- //
// Logic
// -------------------------------------------------------------------------------------------------------------------------- //
// FSM
	always @(posedge clk) begin
		if(~rstn) begin
			c_state <= S_IDLE;
		end
		else begin
			c_state <= n_state;
		end
	end
	always @* begin
		n_state = c_state;
		case (c_state)
			S_IDLE: n_state = S_FETCH;
			S_FETCH: n_state = (cnt == 7)? S_CRC:S_FETCH;
			S_CRC  : n_state = (cnt == 7)? S_END:S_CRC;
			
			//default: 
		endcase
	end

// code info.
	always @(posedge clk) begin
		if(~rstn) begin
			code_reg <= 0;
			mode_reg <= 0;
		end
		else if(set) begin
			code_reg <= code;
			mode_reg <= mode;
		end
	end

// counter
	always @(posedge clk) begin
		if(~rstn) begin
			cnt <= 0;
		end
		else begin
			cnt <= cnt_n;
		end
	end

	always @* begin
		cnt_n = cnt;
		case(c_state)
			S_FETCH: cnt_n = (cnt == 7)? 0:cnt + 1;
		endcase
	end

// llr
	always @(posedge clk) begin
		if(~rstn) begin
			llr <= 0;
		end
		else begin
			llr <= llr_n;
		end
	end
	always @* begin
		llr_n = llr;
		case(c_state)
			S_FETCH: begin
				case(cnt)
					0: llr_n[1:7] = {idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					1: llr_n[8:15] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					2: llr_n[16:23] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					3: llr_n[24:31] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					4: llr_n[32:39] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					5: llr_n[40:47] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					6: llr_n[48:55] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
					7: llr_n[56:63] = {idata[63], idata[55], idata[47], idata[39], idata[31], idata[23], idata[15], idata[7]};
				endcase
			end
		endcase
	end









// output 
	always @* begin
		if(c_state == S_FETCH) ready_reg = 1;
		else ready_reg = 0;
	end

endmodule

